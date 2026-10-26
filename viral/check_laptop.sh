#!/usr/bin/env bash
# check_laptop.sh - can this computer re-run Space Ranger, or should it use the light route (route B)?
# Usage: bash check_laptop.sh [folder where the data/outputs will live]      (changes nothing)
set -uo pipefail
DIR="${1:-.}"

os=$(uname -s); arch=$(uname -m); wsl="no"
if [[ "$os" == "Linux" ]] && grep -qi microsoft /proc/version 2>/dev/null; then wsl="yes (WSL2)"; fi
if [[ "$os" == "Darwin" ]]; then
  cores=$(sysctl -n hw.ncpu)
  ram_gb=$(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))
  avx=$(sysctl -n machdep.cpu.features machdep.cpu.leaf7_features 2>/dev/null | grep -o 'AVX[0-9.]*' | sort -u | tr '\n' ' ')
else
  cores=$(nproc 2>/dev/null || grep -c ^processor /proc/cpuinfo)
  ram_gb=$(awk '/MemTotal/ {printf "%d", $2 / 1024 / 1024}' /proc/meminfo)
  avx=$(grep -o -w -E 'avx2|avx' /proc/cpuinfo 2>/dev/null | sort -u | tr '\n' ' ')
fi
free_gb=$(df -Pk "$DIR" | awk 'NR == 2 {printf "%d", $4 / 1024 / 1024}')
have() { command -v "$1" >/dev/null 2>&1 && echo "yes" || echo "no"; }
pysam=$(python3 -c "import pysam" >/dev/null 2>&1 && echo "yes" || echo "no")

echo "== This computer"
printf "  %-22s %s\n" "OS / architecture" "$os / $arch" "Inside WSL" "$wsl" "CPU cores" "$cores" \
  "RAM" "${ram_gb} GB" "Free disk in $DIR" "${free_gb} GB" "AVX / AVX2" "${avx:-none found}"
printf "  %-22s %s\n" "spaceranger" "$(have spaceranger)" "samtools" "$(have samtools)" \
  "minimap2" "$(have minimap2)" "python3 + pysam" "$pysam"

echo
echo "== Space Ranger re-run (merged tomato + TYLCV reference)"
sr_ok=1
if [[ "$os" == "Darwin" ]]; then
  echo "  NO: Space Ranger runs only on Linux. Use route B on this Mac."; sr_ok=0
elif [[ "$arch" != "x86_64" ]]; then
  echo "  NO: Space Ranger needs an x86_64 (Intel/AMD) CPU; this one is $arch. Use route B."; sr_ok=0
elif [[ -z "$avx" ]]; then
  echo "  NO: the CPU does not report AVX, which Space Ranger requires. Use route B."; sr_ok=0
else
  short=""
  [[ "$cores" -lt 8 ]] && short+=" cores ($cores < 8)"
  [[ "$ram_gb" -lt 64 ]] && short+=" RAM (${ram_gb} GB < 64 GB)"
  [[ "$free_gb" -lt 1000 ]] && short+=" disk (${free_gb} GB < 1 TB)"
  if [[ -n "$short" ]]; then
    echo "  POSSIBLE BUT BELOW 10x's MINIMUM:$short."
    echo "  The tomato genome is smaller than human, so mkref/count may still finish with"
    echo "  --localcores=$cores --localmem=$(( ram_gb * 85 / 100 )) and --no-bam, but this is untested: try ONE capture area first."
  else
    echo "  YES: meets 10x's minimum (8 cores, 64 GB RAM, 1 TB disk). Run one capture area at a time."
  fi
  [[ "$wsl" != "no" ]] && echo "  WSL2: raise its memory cap in %UserProfile%\\.wslconfig (memory=..., processors=...) and keep data inside the Linux filesystem, not /mnt/c."
fi

echo
echo "== Route B (unmapped reads from your original BAMs -> TYLCV)"
missing=""
[[ "$(have samtools)" == "no" ]] && missing+=" samtools"
[[ "$(have minimap2)" == "no" ]] && missing+=" minimap2"
[[ "$pysam" == "no" ]] && missing+=" pysam"
if [[ -n "$missing" ]]; then
  echo "  Install first:$missing"
  if command -v conda >/dev/null 2>&1; then   # bioconda has osx-arm64, osx-64 and linux builds of all three
    echo "  conda create -y -n tylcv_tools -c conda-forge -c bioconda samtools minimap2 pysam && conda activate tylcv_tools"
  elif [[ "$os" == "Darwin" ]]; then          # no conda: Homebrew bottles + pysam wheels
    echo "  brew install samtools minimap2"
    echo "  python3 -m venv ~/tylcv && source ~/tylcv/bin/activate && pip install pysam"
  else
    echo "  conda create -n tylcv -c conda-forge -c bioconda samtools minimap2 pysam python=3.11 && conda activate tylcv"
  fi
else
  echo "  READY: tools found. It needs only a few GB of RAM; disk must hold your BAMs plus small outputs."
fi
[[ "$sr_ok" == 0 ]] && echo "  Recommended on this computer."
