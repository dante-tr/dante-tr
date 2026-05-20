#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

motifs=${1:-}
bam_file=${2:-}
output=${3:-}
script_dir="$(dirname "$(realpath "$0")")"

# %% prerequisities
# compile rust
# cargo build --release
# activate python env
source "$script_dir/dante_py/.venv/bin/activate"
# maybe I should use this instead
# uv run python script.py

# %% run rust part
"$script_dir/target/release/dante_cli"  \
    -m "$motifs"                        \
    -b "$bam_file"                      \
    -o "$output"                        \
    --output-bams

# %% run python part
for rem_output in "$output"/motifs/*.annotations.tsv; do
    if [[ $(wc -l "$rem_output" | cut -f1 -d " ") -eq 1 ]]; then
        rm "$rem_output"
        echo "Oopsie! Zero reads. Skipping..."
        continue
    fi

    py_output=${rem_output%.*}
    mkdir -p "${py_output}"
    python "$script_dir/dante_py/dante_remastr_simple.py"   \
        -i "$rem_output"                                    \
        -o "$py_output"                                     \
        --verbose
done

# %% aggregate python outputs
mkdir -p "${output}/report"
python "$script_dir/dante_py/dante_remastr_simple_agg_and_report.py"    \
    -r "${output}/motifs/"*.annotations.tsv                             \
    -j "${output}/motifs/"*/data.json                                   \
    -o "${output}/report"

