# isolate-genome-assembler-report: Quarto + the Python the report template imports.
# Built by .github/workflows/build-images.yml; locally:
#   docker build -t ghcr.io/timrozday-mgnify/isolate-genome-assembler-report:latest containers/report
FROM python:3.12-slim

ARG QUARTO_VERSION=1.9.37
# Set by buildx; the .deb exists for amd64 and arm64, so the image builds on either.
ARG TARGETARCH

# procps provides `ps`, which Nextflow needs for task metrics.
ADD https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.deb /tmp/quarto.deb
RUN apt-get update \
    && apt-get install -y --no-install-recommends /tmp/quarto.deb procps \
    && rm -rf /var/lib/apt/lists/* /tmp/quarto.deb

# papermill applies `quarto render -P`; ipykernel runs the cells.
RUN pip install --no-cache-dir \
        papermill==2.7.0 \
        ipykernel==7.3.0 \
        pandas==3.0.6 \
        plotly==6.8.0

# Quarto and Jupyter cache under HOME, which is not writable for an arbitrary container
# user; the process also points HOME at its work directory.
ENV HOME=/tmp
