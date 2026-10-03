# TopoQSPR v1.0 command-line image.
#   docker build -t topoqspr .
#   docker run --rm -v "$PWD:/work" -v topoqspr-cache:/root/.cache/topoqspr topoqspr \
#       run --config /work/examples/pilot_migraine.json --out /work/results/pilot
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# git lets the run record name the commit when the source tree is mounted
RUN apt-get update \
    && apt-get install -y --no-install-recommends git libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/topoqspr
COPY requirements-lock.txt ./
RUN pip install -r requirements-lock.txt

COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-deps .

WORKDIR /work
ENTRYPOINT ["topoqspr"]
CMD ["--help"]
