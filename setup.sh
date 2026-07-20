#!/usr/bin/env bash
# REACTR setup script — macOS (Apple Silicon) and Linux (x86_64)
set -euo pipefail

REACTR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OS="$(uname -s)"
ARCH="$(uname -m)"

echo "=== REACTR Setup ==="
echo "Detected: $OS ($ARCH)"
echo ""

# ---------- 1. OS-specific prep ----------
if [[ "$OS" == "Darwin" ]]; then
    echo "--- macOS setup ---"
    if [[ "$ARCH" != "arm64" ]]; then
        echo "Warning: this script targets Apple Silicon (arm64). Intel Macs are untested."
    fi
    xcode-select --install 2>/dev/null || echo "Xcode command line tools already installed."
    softwareupdate --install-rosetta --agree-to-license 2>/dev/null || echo "Rosetta already installed or not applicable."
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"

elif [[ "$OS" == "Linux" ]]; then
    echo "--- Linux setup ---"
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

    # MCScanX needs a C/C++ compiler; not guaranteed present on minimal Linux images/HPC nodes
    if ! command -v g++ &> /dev/null; then
        echo "No system g++ found — will install compilers via conda-forge instead of assuming apt/yum access."
        NEED_CONDA_COMPILERS=1
    else
        echo "System g++ found: $(g++ --version | head -n1)"
        NEED_CONDA_COMPILERS=0
    fi

else
    echo "Unsupported OS: $OS. REACTR currently supports macOS and Linux (see README for Windows via WSL2)."
    exit 1
fi

# ---------- 2. Install Miniforge if conda isn't already available ----------
if ! command -v conda &> /dev/null; then
    echo ""
    echo "Installing Miniforge..."
    curl -L "$MINIFORGE_URL" -o /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
    rm /tmp/miniforge.sh
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
else
    echo "conda already installed, skipping Miniforge install."
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi

# ---------- 3. Create/update the environment ----------
echo ""
echo "Creating conda environment from environment.yaml (this can take a while)..."
if conda env list | grep -q "^reactr "; then
    echo "Environment 'reactr' already exists — updating instead of recreating."
    conda env update --file "$REACTR_DIR/environment.yaml" --prune
else
    conda env create --file "$REACTR_DIR/environment.yaml"
fi
conda activate reactr

# On Linux without a system compiler, pull one in via conda-forge for the MCScanX build step below
if [[ "$OS" == "Linux" && "${NEED_CONDA_COMPILERS:-0}" == "1" ]]; then
    echo "Installing conda-forge compilers into the reactr environment..."
    conda install -y -c conda-forge gxx_linux-64 gcc_linux-64
fi

# ---------- 4. PATH setup ----------
export PATH="$REACTR_DIR/bin:$REACTR_DIR/preset:$PATH"
SHELL_RC="$HOME/.bashrc"
[[ "$SHELL" == *zsh* ]] && SHELL_RC="$HOME/.zshrc"
if ! grep -qs "reactrv1.0/preset" "$SHELL_RC" 2>/dev/null; then
    echo "export PATH=\"$REACTR_DIR/bin:$REACTR_DIR/preset:\$PATH\"" >> "$SHELL_RC"
    echo "Added REACTR bin/ and preset/ to PATH in $SHELL_RC (restart your shell or 'source $SHELL_RC' to pick it up)."
fi

# ---------- 5. Non-conda dependencies ----------
mkdir -p "$REACTR_DIR/preset"
cd "$REACTR_DIR/preset"

# FlashFry
if [[ ! -f "FlashFry.jar" ]]; then
    echo ""
    echo "Downloading FlashFry..."
    wget -q https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar -O FlashFry.jar
else
    echo "FlashFry.jar already present, skipping."
fi

# MCScanX (built from source — needs a C++ compiler, handled above)
if [[ ! -f "$REACTR_DIR/preset/MCScanX" ]]; then
    echo ""
    echo "Building MCScanX..."
    wget -q https://github.com/wyp1125/MCScanX/archive/refs/heads/master.zip -O MCScanX.zip
    unzip -q MCScanX.zip
    (cd MCScanX-master && make)
    mv MCScanX-master/MCScanX MCScanX-master/duplicate_gene_classifier "$REACTR_DIR/preset/"
    rm -rf MCScanX-master MCScanX.zip
else
    echo "MCScanX already built, skipping."
fi

# Pfam-A database
if [[ ! -f "pfam/Pfam-A.hmm" ]]; then
    echo ""
    echo "Downloading Pfam-A database (large file — may take a while)..."
    mkdir -p pfam && cd pfam
    wget -q http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    cd ..
else
    echo "Pfam-A.hmm already present, skipping download."
fi
hmmpress -f pfam/Pfam-A.hmm

cd "$REACTR_DIR"

# ---------- 6. Verify key binaries resolved correctly ----------
echo ""
echo "--- Verifying installation ---"
MISSING=0
for tool in snakemake diamond hmmscan iqtree meme primer3_core blastp gmap; do
    if ! command -v "$tool" &> /dev/null; then
        echo "  WARNING: '$tool' not found on PATH"
        MISSING=1
    fi
done
[[ -f "$REACTR_DIR/preset/MCScanX" ]] || { echo "  WARNING: MCScanX binary not found"; MISSING=1; }
[[ -f "$REACTR_DIR/preset/duplicate_gene_classifier" ]] || { echo "  WARNING: duplicate_gene_classifier binary not found"; MISSING=1; }
[[ -f "$REACTR_DIR/preset/FlashFry.jar" ]] || { echo "  WARNING: FlashFry.jar not found"; MISSING=1; }

if [[ "$MISSING" == "0" ]]; then
    echo "All core dependencies verified."
else
    echo "Some dependencies are missing — check the warnings above before running the pipeline."
fi

# ---------- 7. Wrap-up ----------
echo ""
echo "=== Setup complete ==="
echo "Note: DeepLoc2 (subcellular localization) requires a separate academic-license"
echo "install from https://services.healthtech.dtu.dk/services/DeepLoc-2.0/"
echo "If installed, uncomment the relevant lines in MainPipeline.smk."
echo ""
echo "Next steps:"
echo "  1. conda activate reactr"
echo "  2. Edit config.yaml with your taxon IDs, assembly accessions, and query sequences"
echo "  3. snakemake -s LoadDatasets.smk --cores 8"
echo "  4. snakemake -s MainPipeline.smk --cores all"
