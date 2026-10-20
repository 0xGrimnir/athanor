# Root filesystem of the Athanor Linux engines pack, exported and imported into WSL2.
#   /opt/athanor/venv   Python: PySCF, GPU4PySCF (CUDA 13), RDKit
#   /opt/athanor/bio    bioconda: gargammel, BWA, samtools, AdapterRemoval, mapDamage2
#   /opt/athanor/refs   reference genomes fetched from NCBI at build time (manifest.json)
FROM ubuntu:24.04
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-venv ca-certificates procps curl bzip2 \
 && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /opt/athanor/requirements.txt
RUN python3 -m venv /opt/athanor/venv \
 && /opt/athanor/venv/bin/pip install --no-cache-dir -r /opt/athanor/requirements.txt \
 && /opt/athanor/venv/bin/pip freeze > /opt/athanor/lock.txt \
 && find /opt/athanor/venv -name "__pycache__" -type d -prune -exec rm -rf {} +
COPY bio-environment.yml /opt/athanor/bio-environment.yml
RUN curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local bin/micromamba \
 && micromamba create -y -p /opt/athanor/bio -f /opt/athanor/bio-environment.yml \
 && micromamba list -p /opt/athanor/bio --explicit --md5 > /opt/athanor/bio-lock.txt \
 && micromamba clean -a -y \
 && rm -rf /root/.cache
COPY fetch_refs.py /opt/athanor/fetch_refs.py
RUN /opt/athanor/venv/bin/python /opt/athanor/fetch_refs.py /opt/athanor/refs
COPY wsl.conf /etc/wsl.conf
# The NVIDIA driver's CUDA library is mounted by WSL at /usr/lib/wsl/lib.
ENV LD_LIBRARY_PATH=/usr/lib/wsl/lib
