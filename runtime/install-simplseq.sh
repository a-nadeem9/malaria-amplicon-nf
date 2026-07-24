#!/usr/bin/env bash
set -euo pipefail

VERSION="${SIMPLSEQ_VERSION:-v1.0.3}"
TARBALL="runtime.tar.gz"
CHECKSUMS="SHA256SUMS.txt"
DEFAULT_BASE_URL="https://github.com/a-nadeem9/malaria-amplicon-nf/releases/download/${VERSION}"
BASE_URL="${SIMPLSEQ_INSTALL_BASE_URL:-$DEFAULT_BASE_URL}"
BUNDLED_RUNTIME_DIR="${SIMPLSEQ_BUNDLED_RUNTIME_DIR:-}"

CACHE_DIR="${HOME}/.cache/simplseq/${VERSION}"
SIMPLSEQ_HOME="${HOME}/.local/share/simplseq"
VERSION_DIR="${SIMPLSEQ_HOME}/versions/${VERSION}"
STAGED_VERSION_DIR="${VERSION_DIR}.candidate"
ENV_DIR="${SIMPLSEQ_HOME}/envs/${VERSION}"
LOG_DIR="${SIMPLSEQ_HOME}/logs"
BIN_DIR="${HOME}/.local/bin"
LAUNCHER_PATH="${BIN_DIR}/simplseq"
STAGED_LAUNCHER="${BIN_DIR}/.simplseq-${VERSION}.candidate"
LOG_FILE="${LOG_DIR}/install-${VERSION}.log"
MICROMAMBA="${SIMPLSEQ_HOME}/bin/micromamba"
REUSE_ENV="${SIMPLSEQ_REUSE_ENV:-1}"
DINEMITES_SHA="210e38852a9911d1411f4948f4ab47b46ddd71ed"
PATH_WAS_MISSING=0
UNAME_S="$(uname -s)"
UNAME_M="$(uname -m)"

say() {
  printf '\n== %s ==\n' "$1"
}

banner() {
  cat <<EOF
======================================================
  >_ malaria-amplicon-nf ${VERSION}
     Linux / WSL / macOS browser workflow setup
     Nextflow + Conda/Mamba runtime
======================================================
EOF
}

fail() {
  printf 'ERROR: %s\n' "$1" >&2
  exit 1
}

case "$UNAME_S" in
  Linux)
    PLATFORM_LABEL="Linux / WSL"
    MAMBA_SUBDIR="linux-64"
    CONDA_PLATFORM="${SIMPLSEQ_CONDA_PLATFORM:-linux-64}"
    PROFILE_FILE="${HOME}/.bashrc"
    SHA256_CHECK=(sha256sum -c)
    ;;
  Darwin)
    PLATFORM_LABEL="macOS"
    PROFILE_FILE="${HOME}/.zshrc"
    SHA256_CHECK=(shasum -a 256 -c)
    case "$UNAME_M" in
      arm64)
        MAMBA_SUBDIR="osx-arm64"
        CONDA_PLATFORM="${SIMPLSEQ_CONDA_PLATFORM:-osx-64}"
        ;;
      x86_64)
        MAMBA_SUBDIR="osx-64"
        CONDA_PLATFORM="${SIMPLSEQ_CONDA_PLATFORM:-osx-64}"
        ;;
      *) fail "Unsupported macOS CPU architecture: $UNAME_M" ;;
    esac
    ;;
  *)
    fail "This installer supports Linux/WSL and macOS only."
    ;;
esac

fetch_asset() {
  local name="$1"
  local target="$2"
  if [[ "$BASE_URL" =~ ^https?:// || "$BASE_URL" =~ ^file:// ]]; then
    curl -fsSL --retry 3 --retry-delay 2 --retry-all-errors \
      "${BASE_URL%/}/${name}" -o "$target" \
      || fail "Required release asset '${name}' is unavailable at ${BASE_URL%/}/${name}"
  else
    cp "${BASE_URL%/}/${name}" "$target" \
      || fail "Required local release asset '${name}' is unavailable in ${BASE_URL}"
  fi
}

mkdir -p "$CACHE_DIR" "$SIMPLSEQ_HOME/bin" "$SIMPLSEQ_HOME/versions" "$SIMPLSEQ_HOME/envs" "$LOG_DIR" "$BIN_DIR"
touch "$LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

banner
echo "Platform: $PLATFORM_LABEL ($UNAME_M)"
echo "Micromamba platform: $MAMBA_SUBDIR"
echo "Conda package platform: $CONDA_PLATFORM"
if [[ -n "$BUNDLED_RUNTIME_DIR" ]]; then
  echo "Runtime source: bundled desktop files"
else
  echo "Base URL: $BASE_URL"
fi
echo "Install log: $LOG_FILE"
if [[ "$REUSE_ENV" == "1" ]]; then
  echo "Runtime mode: reuse existing managed environment when present"
else
  echo "Runtime mode: recreate managed environment for a clean install"
fi

if [[ "$UNAME_S" == "Darwin" && "$UNAME_M" == "arm64" && "$CONDA_PLATFORM" == "osx-64" ]]; then
  if ! /usr/bin/arch -x86_64 /usr/bin/true >/dev/null 2>&1; then
    fail "Apple Silicon macOS installs use the Intel conda runtime for DADA2. Install Rosetta first: softwareupdate --install-rosetta --agree-to-license"
  fi
fi

rm -rf "$STAGED_VERSION_DIR"
rm -f "$STAGED_LAUNCHER"
if [[ -n "$BUNDLED_RUNTIME_DIR" ]]; then
  say "Installing bundled app files"
  [[ -f "${BUNDLED_RUNTIME_DIR%/}/main.nf" ]] \
    || fail "Bundled desktop runtime is missing main.nf at ${BUNDLED_RUNTIME_DIR}"
  [[ -f "${BUNDLED_RUNTIME_DIR%/}/environment.yml" ]] \
    || fail "Bundled desktop runtime is missing environment.yml at ${BUNDLED_RUNTIME_DIR}"
  mkdir -p "$STAGED_VERSION_DIR"
  cp -a "${BUNDLED_RUNTIME_DIR%/}/." "$STAGED_VERSION_DIR/"
else
  say "Downloading release files"
  fetch_asset "$TARBALL" "$CACHE_DIR/$TARBALL"
  fetch_asset "$CHECKSUMS" "$CACHE_DIR/$CHECKSUMS"

  say "Verifying checksum"
  tr -d '\r' < "$CACHE_DIR/$CHECKSUMS" > "$CACHE_DIR/${CHECKSUMS}.unix"
  grep "  ${TARBALL}$" "$CACHE_DIR/${CHECKSUMS}.unix" > "$CACHE_DIR/${TARBALL}.sha256" \
    || fail "No checksum entry found for $TARBALL"
  (cd "$CACHE_DIR" && "${SHA256_CHECK[@]}" "${TARBALL}.sha256")

  say "Installing app files"
  TMP_INSTALL="$(mktemp -d)"
  trap 'rm -rf "$TMP_INSTALL"' EXIT
  tar -xzf "$CACHE_DIR/$TARBALL" -C "$TMP_INSTALL"
  EXTRACTED="$(find "$TMP_INSTALL" -mindepth 1 -maxdepth 1 -type d | head -n 1)"
  [[ -n "$EXTRACTED" ]] || fail "Tarball did not contain an app directory."
  cp -a "$EXTRACTED" "$STAGED_VERSION_DIR"
fi

say "Installing micromamba"
if [[ ! -x "$MICROMAMBA" ]]; then
  MM_TMP="$(mktemp -d)"
  curl -fsSL "https://micro.mamba.pm/api/micromamba/${MAMBA_SUBDIR}/latest" -o "$MM_TMP/micromamba.tar.bz2"
  if command -v bzip2 >/dev/null 2>&1; then
    tar -xjf "$MM_TMP/micromamba.tar.bz2" -C "$MM_TMP"
  elif command -v python3 >/dev/null 2>&1; then
    python3 - "$MM_TMP/micromamba.tar.bz2" "$MM_TMP/bin/micromamba" <<'PY'
import pathlib
import shutil
import sys
import tarfile

archive = pathlib.Path(sys.argv[1])
target = pathlib.Path(sys.argv[2])
with tarfile.open(archive, mode="r:bz2") as bundle:
    member = next(
        (item for item in bundle.getmembers() if item.name.rstrip("/").endswith("bin/micromamba")),
        None,
    )
    if member is None or not member.isfile():
        raise SystemExit("Micromamba archive does not contain bin/micromamba")
    source = bundle.extractfile(member)
    if source is None:
        raise SystemExit("Could not read micromamba from its archive")
    target.parent.mkdir(parents=True, exist_ok=True)
    with source, target.open("wb") as output:
        shutil.copyfileobj(source, output)
PY
  else
    fail "Micromamba requires either bzip2 or python3 to unpack its runtime archive."
  fi
  cp "$MM_TMP/bin/micromamba" "$MICROMAMBA"
  chmod +x "$MICROMAMBA"
  rm -rf "$MM_TMP"
fi

say "Creating managed runtime"
export MAMBA_ROOT_PREFIX="${SIMPLSEQ_HOME}/mamba_root"
export CONDA_PKGS_DIRS="${SIMPLSEQ_HOME}/pkgs"
mkdir -p "$MAMBA_ROOT_PREFIX" "$CONDA_PKGS_DIRS"
cd "$STAGED_VERSION_DIR"
if [[ -d "$ENV_DIR" && "$REUSE_ENV" != "1" ]]; then
  echo "Removing existing managed runtime at $ENV_DIR"
  "$MICROMAMBA" remove -y -p "$ENV_DIR" --all || rm -rf "$ENV_DIR"
fi
LOCK_FILE="$STAGED_VERSION_DIR/locks/${CONDA_PLATFORM}-explicit.txt"
if [[ -f "$LOCK_FILE" && "${SIMPLSEQ_USE_LOCK:-1}" != "0" ]]; then
  echo "Using runtime lock: $LOCK_FILE"
  if [[ -x "$ENV_DIR/bin/python" ]]; then
    "$MICROMAMBA" install -y -p "$ENV_DIR" -f "$LOCK_FILE"
  else
    "$MICROMAMBA" create -y -p "$ENV_DIR" -f "$LOCK_FILE"
  fi
else
  echo "No runtime lock selected; resolving from environment.yml"
  if [[ -x "$ENV_DIR/bin/python" ]]; then
    "$MICROMAMBA" install -y --platform "$CONDA_PLATFORM" -p "$ENV_DIR" -f "$STAGED_VERSION_DIR/environment.yml"
  else
    "$MICROMAMBA" create -y --platform "$CONDA_PLATFORM" -p "$ENV_DIR" -f "$STAGED_VERSION_DIR/environment.yml"
  fi
fi

if [[ "$UNAME_S" == "Linux" ]]; then
  for compiler in gcc g++ c++; do
    if [[ ! -x "$ENV_DIR/bin/$compiler" ]]; then
      candidate="$(find "$ENV_DIR/bin" -maxdepth 1 -type f -name "*-${compiler}" | head -n 1)"
      [[ -n "$candidate" ]] || fail "Managed runtime does not contain a ${compiler} compiler."
      ln -s "$(basename "$candidate")" "$ENV_DIR/bin/$compiler"
    fi
  done
  "$ENV_DIR/bin/g++" --version >/dev/null
fi

if [[ -d "$ENV_DIR/bin/cmdstan" ]]; then
  export CMDSTAN="$ENV_DIR/bin/cmdstan"
fi

# Explicit runtime locks can lag behind environment.yml. Excel metadata is an
# advertised input format, so enforce its Python dependency before marking the
# managed runtime ready.
if ! "$ENV_DIR/bin/python" -c 'import openpyxl' >/dev/null 2>&1; then
  say "Installing Excel metadata support"
  "$ENV_DIR/bin/python" -m pip install "openpyxl>=3.1,<4"
fi
"$ENV_DIR/bin/python" -c 'import openpyxl'

"$ENV_DIR/bin/python" -m pip install --no-deps --force-reinstall "$STAGED_VERSION_DIR"
say "Installing downstream R analysis packages"
# Ensure conda compilers (gcc, g++, make) are on PATH for R package builds.
PATH="$ENV_DIR/bin:$PATH" \
R_LIBS_USER="$ENV_DIR/lib/R/library" \
"$ENV_DIR/bin/Rscript" -e '
  required <- c("jsonlite", "dcifer", "instantiate", "patchwork")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    install.packages(
      missing,
      repos = c("https://cloud.r-project.org", "https://mc-stan.org/r-packages/")
    )
  }
'
DINEMITES_MARKER="$ENV_DIR/lib/R/library/dinemites/.simplseq-sha"
DINEMITES_MODEL_DIR="$ENV_DIR/lib/R/library/dinemites/bin/stan"
dinemites_models_ready() {
  [[ -f "$DINEMITES_MARKER" ]] \
    && [[ "$(cat "$DINEMITES_MARKER")" == "$DINEMITES_SHA" ]] \
    && [[ -x "$DINEMITES_MODEL_DIR/model_infection_probabilities_bayesian" ]] \
    && [[ -x "$DINEMITES_MODEL_DIR/model_infection_probabilities_bayesian_drop_out" ]] \
    && [[ -x "$DINEMITES_MODEL_DIR/model_infection_probabilities_clusters" ]]
}
if ! dinemites_models_ready; then
  DINEMITES_TARBALL="$CACHE_DIR/dinemites-${DINEMITES_SHA}.tar.gz"
  curl -fsSL \
    "https://codeload.github.com/WillNickols/dinemites/tar.gz/${DINEMITES_SHA}" \
    -o "$DINEMITES_TARBALL"
  DINEMITES_LDFLAGS="${LDFLAGS:-}"
  if [[ "$UNAME_S" == "Darwin" && "$UNAME_M" == "arm64" && "$CONDA_PLATFORM" == "osx-64" ]]; then
    # Conda's Intel Clang 19 passes a versioned libLTO path that the macOS 26
    # linker rejects. The DINEMITES models do not use LTO, so select the
    # compatible non-plugin linker path when building them under Rosetta.
    DINEMITES_LDFLAGS="${DINEMITES_LDFLAGS:+${DINEMITES_LDFLAGS} }-mlinker-version=0"
  fi
  PATH="$ENV_DIR/bin:$PATH" \
  R_LIBS_USER="$ENV_DIR/lib/R/library" \
  LDFLAGS="$DINEMITES_LDFLAGS" \
    "$ENV_DIR/bin/R" CMD INSTALL \
      --library="$ENV_DIR/lib/R/library" \
      "$DINEMITES_TARBALL"
  printf '%s\n' "$DINEMITES_SHA" > "$DINEMITES_MARKER"
  dinemites_models_ready \
    || fail "DINEMITES installed without all three compiled model executables."
fi

say "Creating launcher"
cat > "$STAGED_LAUNCHER" <<EOF
#!/usr/bin/env bash
set -euo pipefail

SIMPLSEQ_HOME="\${HOME}/.local/share/simplseq"
VERSION="${VERSION}"
PROJECT_ROOT="\${SIMPLSEQ_HOME}/current"
ENV_DIR="\${SIMPLSEQ_HOME}/envs/\${VERSION}"

export SIMPLSEQ_PROJECT_ROOT="\${PROJECT_ROOT}"
export SIMPLSEQ_ENV_DIR="\${ENV_DIR}"
export SIMPLSEQ_VERSION="\${VERSION}"
export CONDA_PREFIX="\${ENV_DIR}"
export PYTHONPATH="\${PROJECT_ROOT}/src\${PYTHONPATH:+:\${PYTHONPATH}}"
export PATH="\${ENV_DIR}/bin:\${PATH}"
[[ -d "\${ENV_DIR}/bin/cmdstan" ]] && export CMDSTAN="\${ENV_DIR}/bin/cmdstan"

# Explicit lock-file installs do not run micromamba's interactive activation
# hook. Source package activation snippets so compiler and CmdStan variables
# match a normal activated environment.
if [[ -d "\${ENV_DIR}/etc/conda/activate.d" ]]; then
  set +u
  for hook in "\${ENV_DIR}"/etc/conda/activate.d/*.sh; do
    [[ -f "\${hook}" ]] && source "\${hook}"
  done
  set -u
fi

exec "\${ENV_DIR}/bin/python" -m simplseq "\$@"
EOF
chmod +x "$STAGED_LAUNCHER"

say "Checking PATH"
if [[ ":$PATH:" != *":$BIN_DIR:"* ]]; then
  PATH_WAS_MISSING=1
  touch "$PROFILE_FILE"
  if ! grep -q 'malaria-amplicon-nf launcher path' "$PROFILE_FILE"; then
    cat >> "$PROFILE_FILE" <<'EOF'

# malaria-amplicon-nf launcher path
case ":$PATH:" in
  *":$HOME/.local/bin:"*) ;;
  *) export PATH="$HOME/.local/bin:$PATH" ;;
esac
EOF
  fi
  echo "$BIN_DIR is not currently on PATH in this shell."
  echo "Open a new shell, or run:"
  echo "  export PATH=\"$BIN_DIR:\$PATH\""
  echo "PATH update written to: $PROFILE_FILE"
fi

say "Verifying malaria-amplicon-nf"
export SIMPLSEQ_PROJECT_ROOT="$STAGED_VERSION_DIR"
export SIMPLSEQ_ENV_DIR="$ENV_DIR"
export SIMPLSEQ_VERSION="$VERSION"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONPATH="$STAGED_VERSION_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
export PATH="$ENV_DIR/bin:$PATH"
"$ENV_DIR/bin/python" -m simplseq --help >/dev/null
"$ENV_DIR/bin/python" -m simplseq check
"$ENV_DIR/bin/python" -m simplseq run-headless --help >/dev/null
"$ENV_DIR/bin/python" -c 'import flask, waitress'

say "Activating verified runtime"
rm -rf "$VERSION_DIR"
mv "$STAGED_VERSION_DIR" "$VERSION_DIR"
printf '%s\n' "$VERSION" > "$VERSION_DIR/.install-ready"
ln -sfn "$VERSION_DIR" "${SIMPLSEQ_HOME}/current"
mv -f "$STAGED_LAUNCHER" "$LAUNCHER_PATH"

say "Setup complete"
if [[ "$PATH_WAS_MISSING" == "1" ]]; then
  cat <<EOF
Start malaria-amplicon-nf now with:

    "$BIN_DIR/simplseq" run

Future terminals can use:

    simplseq run
EOF
else
  cat <<'EOF'
Start malaria-amplicon-nf with:

    simplseq run
EOF
fi
