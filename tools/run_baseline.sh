#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BASELINE_DIR="$PROJECT_ROOT/novae/nova_test/runs/baseline"
LOGDIR="$PROJECT_ROOT/novae/nova_test/logs"

mkdir -p "$LOGDIR"

echo "Building baseline from current ppn template..."
python3 "$SCRIPT_DIR/build_runs.py" --baseline-only

cd "$BASELINE_DIR"
echo "Running baseline..."
START=$(date +%s)
pwd
./ppn.exe | tee log.txt
mv log.txt "$LOGDIR/baseline.log"
END=$(date +%s)
echo "Done"
echo "Elapsed time: $((END-START)) seconds"
