# mdpeptprot: GROMACS + CHARMM36 peptide-protein MD pipeline.
# The image holds code and tools only. Inputs/results live in mounted folders;
# what changes between machines is the job YAML.
FROM nvidia/cuda:12.6.0-base-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      bash gawk curl bzip2 ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*

# micromamba -> one environment with the exact versions validated locally
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
      | tar -xj -C /usr/local/bin --strip-components=1 bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV CONDA_OVERRIDE_CUDA=12.9
RUN micromamba create -y -n md -c conda-forge \
      python=3.11 "gromacs=2026.3=nompi_cuda*" numpy pyyaml pip \
    && /opt/conda/envs/md/bin/pip install --no-cache-dir pdb2pqr==3.8.0 propka==3.5.1 MDAnalysis==2.11.0 \
    && micromamba clean -a -y
ENV PATH=/opt/conda/envs/md/bin:$PATH

# pipeline code (writable by any user, since the container runs as the host user)
WORKDIR /app
COPY scripts/ /app/scripts/
COPY forcefields/ /app/forcefields/
COPY run_md.py /app/run_md.py
RUN chmod -R a+rwX /app
ENV HOME=/tmp

ENTRYPOINT ["python", "/app/run_md.py"]
