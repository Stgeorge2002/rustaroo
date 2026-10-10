//! Build the pangenome from parsed samples.
//!
//! Homology runs first. Synteny then splits paralogs and fusions and adopts
//! fragments, because those are annotation boundaries rather than different
//! genes. The reference for each finished group is chosen only after that, so
//! a fragment cannot outvote the intact copy it was broken off from.

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::cluster::cluster_families;
use crate::error::{Error, Result};
use crate::genetic_code::is_intact;
use crate::model::{Call, Pools};
use crate::reference::{join_fragments, Allele, choose_reference};
use crate::synteny::{
    compress_clusters, freeze_contexts, membership_loci, merge_fragments, refind, split_fusions, split_paralogs,
};
use crate::types::{Gene, GeneMember, MembershipAudit, Pangenome, Params, ParsedSample};

pub fn build(mut samples: Vec<ParsedSample>, params: &Params) -> Result<Pangenome> {
    params.validate()?;
    if samples.is_empty() {
        return Err(Error::new("no samples"));
    }
    samples.sort_by(|a, b| a.name.cmp(&b.name));
    for pair in samples.windows(2) {
        if pair[0].name == pair[1].name {
            return Err(Error::new(format!(
                "two inputs resolved to the sample name '{}'",
                pair[0].name
            )));
        }
    }

    let sample_names: Vec<String> = samples.iter().map(|sample| sample.name.clone()).collect();
    let mut pools = Pools::default();
    let mut calls = Vec::new();
    let mut annotation_ids = Vec::new();
    let mut contigs = Vec::with_capacity(samples.len());
    for (sample_index, sample) in samples.iter().enumerate() {
        contigs.push(
            sample
                .contigs
                .iter()
                .map(|contig| contig.sequence.clone())
                .collect::<Vec<_>>(),
        );
        for gene in &sample.genes {
            if gene.dna.len() < params.min_cds_bp {
                continue;
            }
            let source = if params.audit_members {
                let id = annotation_ids.len() as u32;
                annotation_ids.push(gene.annotation_id.clone());
                Some(id)
            } else {
                None
            };
            calls.push(Call {
                sample: sample_index as u32,
                contig: gene.contig as u32,
                start: gene.start,
                end: gene.end,
                strand: gene.strand,
                dna: pools.intern_dna(gene.dna.clone()),
                protein: pools.intern_protein(gene.protein.clone()),
                cluster: 0,
                context: 0,
                locus: None,
                source,
                fusion_piece: false,
                symbol: gene.symbol.clone(),
                product: gene.product.clone(),
            });
        }
    }
    if calls.is_empty() {
        return Err(Error::new(
            "no coding sequences passed the minimum length",
        ));
    }
    eprintln!(
        "rustroo: {} samples, {} coding sequences, {} unique proteins",
        sample_names.len(),
        calls.len(),
        pools.protein.len()
    );

    // `pools.protein` is already deduplicated by content (interning hashes on
    // insert), so clustering it directly gives one family id per unique
    // protein. Cloning a protein per call first (235k clones on a 167-sample
    // run, for only ~8.6k distinct proteins) made `cluster_families` redo
    // that same deduplication from scratch over the expanded copy.
    eprintln!("rustroo: homology clustering on 1 thread");
    let families = cluster_families(
        &pools.protein,
        params.precluster_identity,
        params.family_identity,
        params.length_ratio,
    );
    for call in calls.iter_mut() {
        call.cluster = families[call.protein as usize];
    }
    compress_clusters(&mut calls);
    let family_count = calls.iter().map(|call| call.cluster).max().map_or(0, |id| id + 1);
    eprintln!("rustroo: {family_count} families after homology");

    eprintln!("rustroo: splitting fusions");
    let fusions = split_fusions(&mut calls, &mut pools);
    eprintln!("rustroo: split {fusions} fusion annotations");
    compress_clusters(&mut calls);
    eprintln!("rustroo: {} groups after fusion repair", dense_family_count(&calls));
    freeze_contexts(&mut calls, &pools);
    let contig_lengths: Vec<Vec<usize>> = contigs.iter()
        .map(|sample| sample.iter().map(String::len).collect()).collect();
    for round in 1..=4 {
        eprintln!("rustroo: synteny round {round}");
        let before = dense_family_count(&calls);
        let split = split_paralogs(&mut calls, &pools);
        let merged = merge_fragments(&mut calls, &pools, &contig_lengths);
        compress_clusters(&mut calls);
        let after = dense_family_count(&calls);
        eprintln!(
            "rustroo: synteny round {round}: {before} -> {} after {split} paralog splits -> {after} after merges",
            before + split
        );
        // Stable anchors and retained loci make splitting idempotent until a
        // real merge changes the evidence; a split-only pass is already settled.
        if !merged {
            eprintln!("rustroo: synteny settled on round {round}");
            break;
        }
        if round == 4 {
            eprintln!("rustroo: synteny reached the 4-round limit while groups were still changing");
        }
    }
    if params.refind {
        eprintln!("rustroo: searching intergenic gaps");
        let found = refind(&mut calls, &mut pools, &contigs);
        eprintln!("rustroo: refound {found} genes in intergenic gaps");
        if found > 0 {
            let split = split_paralogs(&mut calls, &pools);
            compress_clusters(&mut calls);
            eprintln!(
                "rustroo: {split} paralog splits after refinding, {} groups",
                dense_family_count(&calls)
            );
        }
    }

    let mut groups: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, call) in calls.iter().enumerate() {
        groups.entry(call.cluster).or_default().push(index);
    }
    let reps = cluster_reps(&groups, &calls, &pools);
    let locus_ids = membership_loci(&calls, &pools);
    let mut genes = Vec::with_capacity(groups.len());
    let mut ids: Vec<u32> = groups.keys().copied().collect();
    ids.sort_by(|&a, &b| {
        groups[&b]
            .len()
            .cmp(&groups[&a].len())
            .then(a.cmp(&b))
    });
    let total = ids.len();
    eprintln!("rustroo: choosing references for {total} groups");
    // Everything about one group's reference - its alleles, the fusion/
    // extension check, presence and length stats - only reads `calls`,
    // `pools`, and `reps`, none of which are mutated here, so every group is
    // independent of every other and the search can run across threads.
    // Only turning a symbol into a final, de-duplicated gene id has to stay
    // sequential, since a later group's name depends on which names earlier
    // groups (in this exact, already-sorted order) have already taken.
    let progress = AtomicUsize::new(0);
    let partials: Vec<Option<PartialGene>> = ids
        .par_iter()
        .map(|&cluster| {
            let partial = compute_partial_gene(
                cluster,
                &groups[&cluster],
                &calls,
                &pools,
                &reps,
                &locus_ids,
                sample_names.len(),
            );
            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if total >= 100 && (done == 1 || done == total || done % 250 == 0) {
                eprintln!("rustroo: reference {done}/{total}");
            }
            partial
        })
        .collect();
    let mut used_names = HashSet::new();
    let mut name_of_cluster = HashMap::new();
    let mut anonymous = 0u32;
    for (&cluster, partial) in ids.iter().zip(partials) {
        let Some(partial) = partial else { continue };
        let id = unique_name(&partial.symbol, &mut anonymous, &mut used_names);
        if params.audit_members {
            name_of_cluster.insert(cluster, id.clone());
        }
        genes.push(Gene {
            id,
            product: partial.product,
            reference_dna: partial.reference_dna,
            intact: partial.intact,
            selection: partial.selection,
            samples_present: partial.samples_present,
            samples_with_exact_reference: partial.samples_with_exact_reference,
            distinct_alleles: partial.distinct_alleles,
            intact_alleles: partial.intact_alleles,
            presence: partial.presence,
            lengths: partial.lengths,
        });
    }
    genes.sort_by(|a, b| {
        b.samples_present
            .cmp(&a.samples_present)
            .then(a.id.cmp(&b.id))
    });
    let membership_audit = if params.audit_members {
        let gene_indexes: HashMap<&str, usize> = genes.iter().enumerate()
            .map(|(index, gene)| (gene.id.as_str(), index)).collect();
        let mut members = Vec::with_capacity(calls.len());
        for (call, &locus) in calls.iter().zip(&locus_ids) {
            let Some(name) = name_of_cluster.get(&call.cluster) else { continue };
            members.push(GeneMember {
                gene: gene_indexes[name.as_str()],
                sample: call.sample,
                contig: call.contig,
                annotation_id: call.source.map(|id| annotation_ids[id as usize].clone())
                    .unwrap_or_default(),
                start: call.start,
                end: call.end,
                strand: call.strand,
                locus,
                dna: call.dna,
                origin: if call.source.is_none() { "refound" }
                    else if call.fusion_piece { "fusion_piece" } else { "annotated" },
            });
        }
        Some(MembershipAudit {
            contig_names: samples.iter().map(|sample| {
                sample.contigs.iter().map(|contig| contig.name.clone()).collect()
            }).collect(),
            members,
            sequences: pools.dna,
        })
    } else {
        None
    };
    eprintln!("rustroo: {} genes", genes.len());
    Ok(Pangenome {
        sample_names,
        genes,
        membership_audit,
    })
}

fn dense_family_count(calls: &[Call]) -> usize {
    calls.iter().map(|call| call.cluster as usize).max().map_or(0, |id| id + 1)
}

/// Everything about one group's chosen reference except its final, unique
/// gene id - that part alone has to be assigned sequentially afterward, so
/// it is deliberately not computed here. Pulled out of `build` so the search
/// for each group can be an independent unit of parallel work.
struct PartialGene {
    symbol: String,
    product: String,
    reference_dna: String,
    intact: bool,
    selection: &'static str,
    presence: Vec<bool>,
    lengths: Vec<u32>,
    samples_present: u32,
    samples_with_exact_reference: u32,
    distinct_alleles: u32,
    intact_alleles: u32,
}

fn compute_partial_gene(
    cluster: u32,
    indexes: &[usize],
    calls: &[Call],
    pools: &Pools,
    reps: &HashMap<u32, Vec<u8>>,
    locus_ids: &[u32],
    sample_count: usize,
) -> Option<PartialGene> {
    let alleles = alleles_of(indexes, calls, pools);
    let joined = joined_allele(indexes, calls, pools, locus_ids);
    let mut with_join = alleles.clone();
    if let Some(dna) = joined {
        let longer = alleles.iter().all(|allele| dna.len() > allele.dna.len());
        if longer {
            with_join.push(Allele {
                dna,
                count: 0,
                synthetic: true,
            });
        }
    }
    let choice = choose_reference(&with_join, reps, cluster)?;
    let mut presence = vec![false; sample_count];
    let mut lengths = vec![0u32; sample_count];
    let mut intervals: HashMap<(u32, u32), Vec<(u32, u32)>> = HashMap::new();
    let mut exact = HashSet::new();
    for &index in indexes {
        let call = &calls[index];
        presence[call.sample as usize] = true;
        intervals
            .entry((call.sample, call.contig))
            .or_default()
            .push((call.start, call.end));
        if pools.dna[call.dna as usize] == choice.dna {
            exact.insert(call.sample);
        }
    }
    for ((sample, _), spans) in intervals {
        lengths[sample as usize] += covered_bases(&spans);
    }
    let distinct = alleles.len() as u32;
    let intact_alleles = alleles
        .iter()
        .filter(|allele| is_intact(&allele.dna))
        .count() as u32;
    Some(PartialGene {
        symbol: preferred_symbol(indexes, calls),
        product: preferred_product(indexes, calls),
        reference_dna: choice.dna,
        intact: choice.intact,
        selection: choice.selection,
        samples_present: presence.iter().filter(|present| **present).count() as u32,
        samples_with_exact_reference: exact.len() as u32,
        distinct_alleles: distinct,
        intact_alleles,
        presence,
        lengths,
    })
}

fn alleles_of(indexes: &[usize], calls: &[Call], pools: &Pools) -> Vec<Allele> {
    let mut samples_for_dna: HashMap<u32, HashSet<u32>> = HashMap::new();
    for &index in indexes {
        let call = &calls[index];
        samples_for_dna
            .entry(call.dna)
            .or_default()
            .insert(call.sample);
    }
    let mut alleles: Vec<Allele> = samples_for_dna
        .into_iter()
        .map(|(dna, samples)| Allele {
            dna: pools.dna[dna as usize].clone(),
            count: samples.len() as u32,
            synthetic: false,
        })
        .collect();
    alleles.sort_by(|a, b| b.count.cmp(&a.count).then(a.dna.cmp(&b.dna)));
    alleles
}

fn joined_allele(indexes: &[usize], calls: &[Call], pools: &Pools, locus_ids: &[u32]) -> Option<String> {
    let mut by_locus: HashMap<(u32, u32, u32), Vec<usize>> = HashMap::new();
    for &index in indexes {
        let call = &calls[index];
        by_locus.entry((call.sample, call.contig, locus_ids[index])).or_default().push(index);
    }
    let mut joins: HashMap<String, u32> = HashMap::new();
    for group in by_locus.values() {
        if group.len() < 2 {
            continue;
        }
        let mut ordered = group.clone();
        ordered.sort_by_key(|index| calls[*index].start);
        let strand = calls[ordered[0]].strand;
        if ordered.iter().any(|index| calls[*index].strand != strand) {
            continue;
        }
        if strand < 0 {
            ordered.reverse();
        }
        let mut abutting = true;
        for pair in ordered.windows(2) {
            let gap = if strand >= 0 {
                i64::from(calls[pair[1]].start) - i64::from(calls[pair[0]].end) - 1
            } else {
                i64::from(calls[pair[0]].start) - i64::from(calls[pair[1]].end) - 1
            };
            if !(0..=3).contains(&gap) {
                abutting = false;
                break;
            }
        }
        if !abutting {
            continue;
        }
        let parts: Vec<&str> = ordered
            .iter()
            .map(|index| pools.dna[calls[*index].dna as usize].as_str())
            .collect();
        if let Some(joined) = join_fragments(&parts) {
            *joins.entry(joined).or_default() += 1;
        }
    }
    joins
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
        .map(|(dna, _)| dna)
}

fn cluster_reps(
    groups: &HashMap<u32, Vec<usize>>,
    calls: &[Call],
    pools: &Pools,
) -> HashMap<u32, Vec<u8>> {
    let mut reps = HashMap::new();
    for (&cluster, indexes) in groups {
        let protein = indexes
            .iter()
            .map(|index| &pools.protein[calls[*index].protein as usize])
            .max_by_key(|protein| protein.len())
            .cloned()
            .unwrap_or_default();
        reps.insert(cluster, protein);
    }
    reps
}

fn preferred_symbol(indexes: &[usize], calls: &[Call]) -> String {
    let mut votes: HashMap<String, u32> = HashMap::new();
    let mut seen = HashSet::new();
    for &index in indexes {
        let call = &calls[index];
        if !seen.insert((call.sample, call.symbol.clone())) {
            continue;
        }
        if usable_symbol(&call.symbol) {
            *votes.entry(call.symbol.clone()).or_default() += 1;
        }
    }
    votes
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
        .map(|(symbol, _)| symbol)
        .unwrap_or_default()
}

fn preferred_product(indexes: &[usize], calls: &[Call]) -> String {
    let mut votes: HashMap<String, u32> = HashMap::new();
    for &index in indexes {
        let product = calls[index].product.trim();
        if product.is_empty() {
            continue;
        }
        *votes.entry(product.to_string()).or_default() += 1;
    }
    votes
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
        .map(|(product, _)| product)
        .unwrap_or_else(|| "hypothetical protein".to_string())
}

fn usable_symbol(symbol: &str) -> bool {
    let trimmed = symbol.trim();
    if trimmed.is_empty() {
        return false;
    }
    let lower = trimmed.to_ascii_lowercase();
    !matches!(
        lower.as_str(),
        "hypothetical" | "hypothetical protein" | "unknown" | "uncharacterized protein"
    )
}

fn unique_name(symbol: &str, anonymous: &mut u32, used: &mut HashSet<String>) -> String {
    let mut base = sanitize(symbol);
    if !usable_symbol(&base) {
        *anonymous += 1;
        base = format!("group_{anonymous:05}");
    }
    let mut name = base.clone();
    let mut suffix = 2u32;
    while !used.insert(name.clone()) {
        name = format!("{base}_{suffix}");
        suffix += 1;
    }
    name
}

fn sanitize(symbol: &str) -> String {
    let mut out = String::with_capacity(symbol.len());
    for character in symbol.chars() {
        if character.is_ascii_alphanumeric() || matches!(character, '.' | '_' | '-') {
            out.push(character);
        } else {
            out.push('_');
        }
    }
    out
}

fn covered_bases(intervals: &[(u32, u32)]) -> u32 {
    let mut spans = intervals.to_vec();
    spans.sort_by_key(|span| span.0);
    let mut total = 0u32;
    let mut current: Option<(u32, u32)> = None;
    for (start, end) in spans {
        if end < start {
            continue;
        }
        if let Some((open, close)) = current {
            if start <= close.saturating_add(1) {
                current = Some((open, close.max(end)));
            } else {
                total = total.saturating_add(close - open + 1);
                current = Some((start, end));
            }
        } else {
            current = Some((start, end));
        }
    }
    if let Some((open, close)) = current {
        total = total.saturating_add(close - open + 1);
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::{back_translate, translate};
    use crate::types::{Contig, ParsedGene, ParsedSample};

    fn diverse(length: usize, salt: usize) -> Vec<u8> {
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let mut protein = Vec::with_capacity(length);
        protein.push(b'M');
        for index in 1..length {
            let mut residue = alphabet[(index * 5 + salt) % alphabet.len()];
            if residue == b'M' {
                residue = b'L';
            }
            protein.push(residue);
        }
        protein
    }

    fn random_protein(length: usize, mut state: u64) -> Vec<u8> {
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let mut sequence = vec![b'M'];
        for _ in 1..length {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let aa = alphabet[(state % 20) as usize];
            sequence.push(if aa == b'M' { b'L' } else { aa });
        }
        sequence
    }

    fn coding(protein: &[u8], contig: usize, start: u32, symbol: &str) -> ParsedGene {
        let dna = back_translate(protein);
        let end = start + dna.len() as u32 - 1;
        ParsedGene {
            annotation_id: format!("cds_{contig}_{start}"),
            contig,
            start,
            end,
            strand: 1,
            protein: translate(&dna),
            dna,
            symbol: symbol.to_string(),
            product: symbol.to_string(),
        }
    }

    fn sample(name: &str, genes: Vec<ParsedGene>) -> ParsedSample {
        let mut end = 1u32;
        for gene in &genes {
            end = end.max(gene.end + 10);
        }
        ParsedSample {
            name: name.to_string(),
            contigs: vec![Contig {
                name: "contig".to_string(),
                sequence: "N".repeat(end as usize + 10),
            }],
            genes,
        }
    }

    fn params() -> Params {
        Params {
            min_cds_bp: 30,
            refind: false,
            ..Params::default()
        }
    }

    #[test]
    fn membership_audit_does_not_change_references_or_grouping() {
        let protein = diverse(40, 4);
        let other = diverse(40, 7);
        let inputs = vec![
            sample("b", vec![coding(&protein, 0, 300, "gene")]),
            sample("a", vec![coding(&other, 0, 100, "extra"), coding(&protein, 0, 400, "gene")]),
        ];
        let normal = build(inputs.clone(), &params()).unwrap();
        assert!(normal.membership_audit.is_none());
        let mut audit_params = params();
        audit_params.audit_members = true;
        let audited = build(inputs, &audit_params).unwrap();
        assert_eq!(normal.sample_names, audited.sample_names);
        assert_eq!(normal.genes.len(), audited.genes.len());
        for (a, b) in normal.genes.iter().zip(&audited.genes) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.reference_dna, b.reference_dna);
            assert_eq!(a.presence, b.presence);
            assert_eq!(a.lengths, b.lengths);
            assert_eq!(a.selection, b.selection);
        }
        let audit = audited.membership_audit.as_ref().unwrap();
        assert_eq!(audit.members.len(), 3);
        for member in &audit.members {
            let extra = member.start == 100;
            assert_eq!(audited.genes[member.gene].id, if extra { "extra" } else { "gene" });
            let expected = if extra { 100 }
                else if audited.sample_names[member.sample as usize] == "a" { 400 } else { 300 };
            assert_eq!(member.start, expected);
            assert_eq!(member.annotation_id, format!("cds_0_{expected}"));
            assert_eq!(member.origin, "annotated");
            assert_eq!(audit.sequences[member.dna as usize], back_translate(if extra { &other } else { &protein }));
        }
        assert_ne!(audit.members[0].locus, audit.members[1].locus);
    }

    #[test]
    fn rare_full_gene_is_the_reference_and_truncations_count_as_present() {
        let full = diverse(70, 1);
        let short = full[..30].to_vec();
        let flank_l = diverse(40, 8);
        let flank_r = diverse(40, 15);
        let mut samples = Vec::new();
        samples.push(sample(
            "intact",
            vec![
                coding(&flank_l, 0, 1, "flkL"),
                coding(&full, 0, 400, "coreA"),
                coding(&flank_r, 0, 900, "flkR"),
            ],
        ));
        for index in 0..4 {
            samples.push(sample(
                &format!("trunc{index}"),
                vec![
                    coding(&flank_l, 0, 1, "flkL"),
                    coding(&short, 0, 400, "coreA"),
                    coding(&flank_r, 0, 900, "flkR"),
                ],
            ));
        }
        let built = build(samples, &params()).unwrap();
        let gene = built
            .genes
            .iter()
            .find(|gene| gene.id == "coreA")
            .expect("coreA");
        assert_eq!(gene.reference_dna, back_translate(&full));
        assert!(gene.intact);
        assert!(gene.presence.iter().all(|present| *present));
        assert_eq!(gene.samples_with_exact_reference, 1);
    }

    #[test]
    fn contig_edge_alleles_keep_the_full_reference_and_every_samples_membership() {
        let full = random_protein(120, 42);
        let upstream = random_protein(36, 101);
        let downstream = random_protein(36, 202);
        let mut inputs = Vec::new();
        let mut expected_lengths = Vec::new();
        for index in 0..10 {
            let reverse = index % 2 == 1;
            let complete = index < 2;
            let protein = if complete { full.as_slice() } else { &full[30..] };
            let mut middle = coding(protein, 0, if !complete && !reverse { 19 } else { 400 }, "gene");
            if !complete {
                middle.dna.replace_range(0..3, "ATG");
                middle.protein = translate(&middle.dna);
            }
            middle.strand = if reverse { -1 } else { 1 };
            expected_lengths.push(middle.dna.len() as u32);
            let mut genes = Vec::new();
            if complete {
                genes.push(coding(if reverse { &downstream } else { &upstream }, 0, 100, "flank"));
                genes.push(coding(if reverse { &upstream } else { &downstream }, 0, 1100, "flank"));
            } else if reverse {
                genes.push(coding(&downstream, 0, 100, "flank"));
            } else {
                genes.push(coding(&downstream, 0, middle.end + 80, "flank"));
            }
            genes.push(middle);
            inputs.push(sample(&format!("s{index:02}"), genes));
        }
        let mut params = params();
        params.audit_members = true;
        let built = build(inputs, &params).unwrap();
        assert_eq!(built.genes.len(), 3);
        let gene = built.genes.iter().find(|gene| gene.id == "gene").unwrap();
        assert!(gene.intact);
        assert_eq!(gene.reference_dna, back_translate(&full), "common truncations must not define the reference");
        assert_eq!(gene.presence, vec![true; 10]);
        assert_eq!(gene.lengths, expected_lengths);
        assert_eq!(gene.samples_with_exact_reference, 2);
        let audit = built.membership_audit.as_ref().unwrap();
        let members: Vec<_> = audit.members.iter().filter(|member| built.genes[member.gene].id == "gene").collect();
        assert_eq!(members.len(), 10);
        assert!(members.iter().all(|member| !member.annotation_id.is_empty() && member.origin == "annotated"));
        assert_eq!(members.iter().map(|member| member.locus).collect::<HashSet<_>>().len(), 10);
    }

    fn pooled_test_call(
        pools: &mut Pools, protein: &[u8], contig: u32, start: u32, strand: i8,
    ) -> Call {
        let dna = back_translate(protein);
        Call {
            sample: 0, contig, start, end: start + dna.len() as u32 - 1, strand,
            protein: pools.intern_protein(translate(&dna)), dna: pools.intern_dna(dna),
            cluster: 0, context: 0, locus: None, source: Some(0), fusion_piece: false,
            symbol: String::new(), product: String::new(),
        }
    }

    #[test]
    fn reference_joins_use_one_supported_locus_on_either_strand() {
        let full = random_protein(80, 42);
        for reverse in [false, true] {
            let mut pools = Pools::default();
            let strand = if reverse { -1 } else { 1 };
            let first = if reverse { &full[40..] } else { &full[..40] };
            let second = if reverse { &full[..40] } else { &full[40..] };
            let a = pooled_test_call(&mut pools, first, 0, 100, strand);
            let b = pooled_test_call(&mut pools, second, 0, a.end + 1, strand);
            let calls = vec![a, b];
            let loci = membership_loci(&calls, &pools);
            assert_eq!(loci[0], loci[1]);
            assert_eq!(joined_allele(&[0, 1], &calls, &pools, &loci), Some(back_translate(&full)));
        }
    }

    #[test]
    fn reference_joins_reject_cross_contig_copies_full_tandems_and_overlaps() {
        let full = random_protein(80, 42);
        for other_contig in [false, true] {
            let mut pools = Pools::default();
            let a = pooled_test_call(&mut pools, &full[..40], 0, 100, 1);
            let b = pooled_test_call(&mut pools, &full[40..], u32::from(other_contig), a.end + 1, 1);
            let calls = vec![a, b];
            if other_contig {
                let loci = membership_loci(&calls, &pools);
                assert_ne!(loci[0], loci[1]);
                assert!(joined_allele(&[0, 1], &calls, &pools, &loci).is_none());
            } else {
                assert!(joined_allele(&[0, 1], &calls, &pools, &[0, 1]).is_none(),
                    "abutting coordinates alone cannot license a join");
            }
        }
        for overlap in [false, true] {
            let mut pools = Pools::default();
            let a = pooled_test_call(&mut pools, &full, 0, 100, 1);
            let b = pooled_test_call(&mut pools, &full, 0, if overlap { 101 } else { a.end + 1 }, 1);
            let calls = vec![a, b];
            let loci = membership_loci(&calls, &pools);
            assert!(joined_allele(&[0, 1], &calls, &pools, &loci).is_none());
        }
    }

    #[test]
    fn unresolved_copies_on_different_contigs_keep_their_covered_lengths() {
        let protein = random_protein(40, 42);
        let first = coding(&protein, 0, 100, "dup");
        let second = coding(&protein, 1, 100, "dup");
        let bp = first.dna.len() as u32;
        let input = ParsedSample {
            name: "one".to_string(),
            contigs: (0..2).map(|index| Contig {
                name: format!("contig{index}"), sequence: "N".repeat(1000),
            }).collect(),
            genes: vec![first, second],
        };
        let mut params = params();
        params.audit_members = true;
        let built = build(vec![input], &params).unwrap();
        assert_eq!(built.genes.len(), 1, "indistinguishable copies stay explicitly unresolved");
        assert_eq!(built.genes[0].reference_dna, back_translate(&protein));
        assert_eq!(built.genes[0].lengths, vec![bp * 2]);
        let members = &built.membership_audit.as_ref().unwrap().members;
        assert_eq!(members.len(), 2);
        assert_ne!(members[0].locus, members[1].locus);
    }

    #[test]
    fn moderate_terminal_alleles_preserve_full_reference_and_sample_lengths() {
        let full = diverse(100, 1);
        let short = full[..90].to_vec();
        let flank_l = diverse(36, 8);
        let flank_r = diverse(36, 15);
        let make_sample = |name: &str, protein: &[u8]| {
            sample(name, vec![
                coding(&flank_l, 0, 1, "flkL"),
                coding(protein, 0, 400, "coreA"),
                coding(&flank_r, 0, 900, "flkR"),
            ])
        };
        let built = build(
            vec![make_sample("intact", &full), make_sample("short", &short)],
            &params(),
        ).unwrap();
        assert_eq!(built.genes.len(), 3);
        let gene = built.genes.iter().find(|gene| gene.id == "coreA").expect("coreA");
        assert_eq!(gene.reference_dna, back_translate(&full));
        assert!(gene.intact);
        assert_eq!(gene.presence, vec![true, true]);
        assert_eq!(gene.samples_with_exact_reference, 1);
        assert_eq!(gene.lengths, vec![
            back_translate(&full).len() as u32, back_translate(&short).len() as u32,
        ]);
    }

    #[test]
    fn fragments_are_the_same_gene_as_the_intact_copy() {
        let full = diverse(80, 3);
        let left = full[..28].to_vec();
        let right = full[28..].to_vec();
        let flank_l = diverse(36, 8);
        let flank_r = diverse(36, 16);
        let intact = sample(
            "intact",
            vec![
                coding(&flank_l, 0, 1, "flkL"),
                coding(&full, 0, 500, "geneA"),
                coding(&flank_r, 0, 1100, "flkR"),
            ],
        );
        let broken_left = coding(&left, 0, 500, "geneA");
        let mut broken_right = coding(&right, 0, broken_left.end + 1, "geneA");
        broken_right.start = broken_left.end + 1;
        broken_right.end = broken_right.start + broken_right.dna.len() as u32 - 1;
        let broken = sample(
            "broken",
            vec![
                coding(&flank_l, 0, 1, "flkL"),
                broken_left,
                broken_right,
                coding(&flank_r, 0, 1100, "flkR"),
            ],
        );
        let built = build(vec![intact, broken], &params()).unwrap();
        let gene = built.genes.iter().find(|gene| gene.id == "geneA").unwrap();
        assert_eq!(gene.reference_dna, back_translate(&full));
        assert_eq!(gene.samples_present, 2);
    }

    #[test]
    fn paralogs_become_two_genes() {
        let paralog = diverse(40, 4);
        let left = diverse(36, 6);
        let mid = diverse(36, 9);
        let right = diverse(36, 12);
        let copy_a = coding(&paralog, 0, 300, "dup");
        let middle = coding(&mid, 0, copy_a.end + 200, "mid");
        let copy_b = coding(&paralog, 0, middle.end + 200, "dup");
        let built = build(
            vec![sample(
                "one",
                vec![
                    coding(&left, 0, 1, "left"),
                    copy_a,
                    middle,
                    copy_b,
                    coding(&right, 0, 2000, "right"),
                ],
            )],
            &params(),
        )
        .unwrap();
        let copies: Vec<_> = built
            .genes
            .iter()
            .filter(|gene| gene.reference_dna == back_translate(&paralog))
            .collect();
        assert_eq!(copies.len(), 2, "paralogs stayed in one group");
        assert!(copies.iter().all(|gene| gene.samples_present == 1));
    }

    #[test]
    fn touching_full_length_paralogs_are_not_joined_as_fragments() {
        let protein = diverse(40, 4);
        let first = coding(&protein, 0, 300, "dup");
        let second = coding(&protein, 0, first.end + 1, "dup");
        let built = build(vec![sample("one", vec![
            coding(&diverse(36, 6), 0, 1, "left"),
            first,
            second,
            coding(&diverse(36, 12), 0, 1000, "right"),
        ])], &params()).unwrap();
        let copies: Vec<_> = built.genes.iter()
            .filter(|gene| gene.reference_dna == back_translate(&protein)).collect();
        assert_eq!(copies.len(), 2);
        assert!(copies.iter().all(|gene| gene.samples_present == 1));
        assert!(copies.iter().all(|gene| gene.lengths == vec![back_translate(&protein).len() as u32]));
    }

    #[test]
    fn fragmented_paralog_preserves_two_full_references_and_presence() {
        let full = random_protein(80, 42);
        let left_a = random_protein(36, 101);
        let right_a = random_protein(36, 202);
        let left_b = random_protein(36, 303);
        let right_b = random_protein(36, 404);
        let mut samples = Vec::new();
        for index in 0..3 {
            let mut genes = vec![coding(&left_a, 0, 1, "leftA")];
            if index == 1 {
                genes.push(coding(&full, 0, 400, "dupA"));
            } else {
                let head = coding(&full[..40], 0, 400, "dupA");
                let tail = coding(&full[40..], 0, head.end + 1, "dupA");
                genes.extend([head, tail]);
            }
            genes.extend([
                coding(&right_a, 0, 900, "rightA"),
                coding(&left_b, 0, 1100, "leftB"),
                coding(&full, 0, 1400, "dupB"),
                coding(&right_b, 0, 2000, "rightB"),
            ]);
            samples.push(sample(&format!("sample{index}"), genes));
        }
        let built = build(samples, &params()).unwrap();
        assert_eq!(built.genes.len(), 6);
        let copies: Vec<_> = built.genes.iter()
            .filter(|gene| gene.reference_dna == back_translate(&full)).collect();
        assert_eq!(copies.len(), 2);
        assert!(copies.iter().all(|gene| gene.intact && gene.presence == vec![true; 3]));
        let lengths: Vec<_> = copies.iter().map(|gene| gene.lengths.clone()).collect();
        let full_bp = back_translate(&full).len() as u32;
        assert!(lengths.contains(&vec![full_bp; 3]));
        assert!(lengths.contains(&vec![full_bp + 3, full_bp, full_bp + 3]));
    }

    #[test]
    fn a_fusion_is_recorded_as_both_genes() {
        let gene_a = diverse(40, 2);
        let gene_b = diverse(36, 14);
        let mut fusion = gene_a.clone();
        fusion.extend_from_slice(&gene_b[1..]);
        let mut samples = Vec::new();
        for index in 0..3 {
            samples.push(sample(
                &format!("a{index}"),
                vec![
                    coding(&diverse(30, 5), 0, 1, "flk"),
                    coding(&gene_a, 0, 400, "aaa"),
                ],
            ));
            samples.push(sample(
                &format!("b{index}"),
                vec![coding(&gene_b, 0, 400, "bbb")],
            ));
        }
        samples.push(sample(
            "fused",
            vec![
                coding(&diverse(30, 5), 0, 1, "flk"),
                coding(&fusion, 0, 400, "aaa"),
            ],
        ));
        let mut audit_params = params();
        audit_params.audit_members = true;
        let built = build(samples, &audit_params).unwrap();
        let aaa = built.genes.iter().find(|gene| gene.id == "aaa").unwrap();
        let bbb = built.genes.iter().find(|gene| gene.id == "bbb").unwrap();
        assert_eq!(aaa.reference_dna, back_translate(&gene_a));
        assert_eq!(bbb.reference_dna, back_translate(&gene_b));
        let fused = built
            .sample_names
            .iter()
            .position(|name| name == "fused")
            .unwrap();
        assert!(aaa.presence[fused]);
        assert!(bbb.presence[fused]);
        assert!(built
            .genes
            .iter()
            .all(|gene| gene.reference_dna != back_translate(&fusion)));
        let fusion_members: Vec<_> = built.membership_audit.as_ref().unwrap().members.iter()
            .filter(|member| member.sample as usize == fused && member.origin == "fusion_piece")
            .collect();
        assert_eq!(fusion_members.len(), 2);
        assert!(fusion_members.iter().all(|member| member.annotation_id == "cds_0_400"));
        assert_ne!(fusion_members[0].gene, fusion_members[1].gene);
    }

    #[test]
    fn a_missed_gene_is_refound_between_the_same_flanks() {
        let flank_l = diverse(36, 7);
        let gene = diverse(40, 3);
        let flank_r = diverse(36, 18);
        let gene_dna = back_translate(&gene);
        let annotated = |name: &str| {
            let left = coding(&flank_l, 0, 1, "flkL");
            let middle = coding(&gene, 0, left.end + 80, "missed");
            let right = coding(&flank_r, 0, middle.end + 80, "flkR");
            let mut sequence = vec![b'N'; right.end as usize + 10];
            paint(&mut sequence, &left);
            paint(&mut sequence, &middle);
            paint(&mut sequence, &right);
            sample_with(name, sequence, vec![left, middle, right])
        };
        let left = coding(&flank_l, 0, 1, "flkL");
        // The CDS sits inside the gap, unannotated. The flanks are the same
        // pair that brackets it in the other samples, which is the only
        // evidence refind is allowed to use.
        let gene_start = left.end + 41;
        let right = coding(
            &flank_r,
            0,
            gene_start + gene_dna.len() as u32 + 48,
            "flkR",
        );
        let mut sequence = vec![b'N'; right.end as usize + 10];
        paint(&mut sequence, &left);
        let at = (gene_start - 1) as usize;
        sequence[at..at + gene_dna.len()].copy_from_slice(gene_dna.as_bytes());
        paint(&mut sequence, &right);
        let gap = right.start as i32 - left.end as i32 - 1;
        assert!((60..=12_000).contains(&gap), "{gap}");

        let mut params = params();
        params.refind = true;
        params.audit_members = true;
        let built = build(
            vec![
                annotated("keepA"),
                annotated("keepB"),
                sample_with("gap", sequence, vec![left, right]),
            ],
            &params,
        )
        .unwrap();
        let row = built
            .genes
            .iter()
            .find(|gene| gene.id == "missed")
            .expect("missed");
        let gap_sample = built
            .sample_names
            .iter()
            .position(|name| name == "gap")
            .unwrap();
        assert!(row.presence[gap_sample]);
        assert_eq!(row.lengths[gap_sample], gene_dna.len() as u32);
        assert_eq!(row.reference_dna, gene_dna);
        assert_eq!(row.samples_present, 3);
        let found: Vec<_> = built.membership_audit.as_ref().unwrap().members.iter()
            .filter(|member| member.origin == "refound").collect();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].sample as usize, gap_sample);
        assert!(found[0].annotation_id.is_empty());
        assert_eq!(built.genes[found[0].gene].id, "missed");
    }

    fn paint(sequence: &mut [u8], gene: &ParsedGene) {
        let start = (gene.start - 1) as usize;
        let dna = gene.dna.as_bytes();
        sequence[start..start + dna.len()].copy_from_slice(dna);
    }

    fn sample_with(name: &str, sequence: Vec<u8>, genes: Vec<ParsedGene>) -> ParsedSample {
        ParsedSample {
            name: name.to_string(),
            contigs: vec![Contig {
                name: "contig".to_string(),
                sequence: String::from_utf8(sequence).expect("contig is ascii"),
            }],
            genes,
        }
    }
}
