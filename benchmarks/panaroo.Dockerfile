ARG PANAROO_BASE=staphb/panaroo:1.5.1
FROM ${PANAROO_BASE}

USER root
RUN apt-get update \
    && apt-get install --no-install-recommends --yes time \
    && apt-get clean \
    && /usr/bin/time --version
