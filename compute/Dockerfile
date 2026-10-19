# Root filesystem of the Athanor compute pack, exported and imported into WSL2.
FROM ubuntu:24.04
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-venv ca-certificates \
 && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /opt/athanor/requirements.txt
RUN python3 -m venv /opt/athanor/venv \
 && /opt/athanor/venv/bin/pip install --no-cache-dir -r /opt/athanor/requirements.txt \
 && /opt/athanor/venv/bin/pip freeze > /opt/athanor/lock.txt \
 && find /opt/athanor/venv -name "__pycache__" -type d -prune -exec rm -rf {} +
COPY wsl.conf /etc/wsl.conf
# The NVIDIA driver's CUDA library is mounted by WSL at /usr/lib/wsl/lib.
ENV LD_LIBRARY_PATH=/usr/lib/wsl/lib
