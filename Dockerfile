# TopoQSPR v1.0 command-line image.
#   docker build -t topoqspr \
#       --build-arg GIT_COMMIT="$(git rev-parse HEAD)" \
#       --build-arg GIT_DIRTY="$(git status --porcelain --untracked-files=no | grep -q . && echo 1 || echo 0)" .
#   docker run --rm -v "$PWD:/work" -v topoqspr-cache:/root/.cache/topoqspr topoqspr \
#       run --config /work/examples/pilot_migraine.json --out /work/results/pilot
#
# The package is installed into site-packages, outside any git tree, so the run
# record cannot read the commit from git at run time (mounting the source tree
# does not change that). The build therefore requires GIT_COMMIT and GIT_DIRTY
# and writes them into topoqspr/_build_info.py; provenance.git_info() falls back
# to that file and the record shows code.source = "build_info".
FROM python:3.12-slim

ARG GIT_COMMIT=""
ARG GIT_DIRTY=""

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/topoqspr
# the lock pins every runtime dependency of pyproject.toml (incl. joblib and
# openpyxl for expdata's mp_bradley workbook); the package goes in with --no-deps
COPY requirements-lock.txt ./
RUN pip install -r requirements-lock.txt

COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN echo "$GIT_COMMIT" | grep -Eq '^[0-9a-f]{40}$' \
        || { echo 'build with --build-arg GIT_COMMIT="$(git rev-parse HEAD)" (40 hex digits)' >&2; exit 1; } \
    && case "$GIT_DIRTY" in 0|1) ;; *) echo 'build with --build-arg GIT_DIRTY=0 or 1' >&2; exit 1;; esac \
    && printf '"""Written by the Dockerfile at build time; see provenance.build_info."""\nGIT_COMMIT = "%s"\nGIT_DIRTY = "%s"\n' \
        "$GIT_COMMIT" "$GIT_DIRTY" > src/topoqspr/_build_info.py \
    && pip install --no-deps .
LABEL org.opencontainers.image.revision="$GIT_COMMIT"

WORKDIR /work
ENTRYPOINT ["topoqspr"]
CMD ["--help"]
