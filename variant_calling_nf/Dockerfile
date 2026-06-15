# Container image for variant_calling_nf, built straight from the pinned conda
# env. Base: micromamba (lightweight mamba) — the env in
# envs/variant_calling_env.yml is solved at build time so the image matches
# `-profile conda` exactly. This bakes in the DNA alignment + GATK toolchain
# (bwa, samtools, gatk4, bcftools, fastp, seqtk) plus the Python stack, so the
# germline-calling DAG runs the same anywhere.
#
# Build:  docker build -t variant_calling_nf:0.1.0 .
# Run  :  nextflow run main.nf -profile docker
FROM mambaorg/micromamba:1.5.8

# Solve the pinned environment into the image's base env.
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/variant_calling_env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml \
    && micromamba clean --all --yes

# Make the env's tools the default for every RUN/CMD and for Nextflow's process
# commands (Nextflow invokes a bash login shell inside the container).
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PYTHONNOUSERSITE=1

# Pipeline scripts live in bin/; copy them so the image is self-contained.
COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ /opt/pipeline/bin/
ENV PATH=/opt/pipeline/bin:$PATH
