#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test "$(uname -m)" = x86_64 || {
    echo "MutCode requires a Linux x86_64 build host." >&2
    exit 1
}
command -v apptainer >/dev/null
OUTPUT="${1:-$ROOT/out/containers}"
mkdir -p "$OUTPUT"
OUTPUT="$(cd -- "$OUTPUT" && pwd)"
SNAPSHOT="$(mktemp -d "${TMPDIR:-/tmp}/mutcode-build.XXXXXXXX")"
printf 'Source snapshot: %s\n' "$SNAPSHOT"

cd "$ROOT"
git ls-files --cached --others --exclude-standard -z -- \
    Cargo.toml Cargo.lock src tests containers/mutcode.def | \
    tar --null -T - -cf - | tar -xf - -C "$SNAPSHOT"

cd "$SNAPSHOT"
test -f Cargo.toml
test -f Cargo.lock
test -f containers/mutcode.def
{
    sha256sum Cargo.toml Cargo.lock containers/mutcode.def
    find src tests -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum
} > source.sha256
FINGERPRINT="$(sha256sum source.sha256 | cut -d ' ' -f 1)"
IMAGE="$OUTPUT/mutcode-${FINGERPRINT:0:12}-linux-amd64.sif"
if [ -e "$IMAGE" ]; then
    echo "Image already exists; refusing to overwrite: $IMAGE" >&2
    exit 1
fi
printf 'Source fingerprint: %s\nImage: %s\n' "$FINGERPRINT" "$IMAGE"
cp containers/mutcode.def "$IMAGE.def"
apptainer build --fakeroot \
    --build-arg "SOURCE_SHA256=$FINGERPRINT" \
    "$IMAGE" containers/mutcode.def 2>&1 | tee "$IMAGE.build.log"
apptainer test "$IMAGE"
apptainer exec --cleanenv "$IMAGE" mutcode --version
cd "$OUTPUT"
sha256sum "$(basename "$IMAGE")" > "$(basename "$IMAGE").sha256"
printf 'Built image: %s\nRetained source snapshot: %s\n' "$IMAGE" "$SNAPSHOT"
