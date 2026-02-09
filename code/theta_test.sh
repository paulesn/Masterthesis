#!/usr/bin/env bash
set -euo pipefail

# Change these if your binaries live elsewhere
A_BIN="${A_BIN:-./src/daniel/theta-graph/cmake-build-release/createThetaGraph}"   # C++ program A
B_BIN="${B_BIN:-./cmake-build-release/analyse}"   # C++ program B
#B_BIN="${B_BIN:-./cmake-build-release/analyseVis}"
#base_graph="./../data/Archipelago.map-coastlines.txt.graph"
base_graph="./../data/0025.32.fmi"


# The sequence of i values
start=24
end=25
k_list=(24 32 50 64 75 80 90 100 124 128)
# ls "../data"


#for i in "${k_list[@]}"; do
for ((i=start; i<=end; i++)); do
  echo "Running Daniel's Theta spanner generation with theta=$i..."
  "$A_BIN" \
    -g $base_graph \
    -o "./../data/daniel-25.fmi" \
    -n "$i"

  echo "A finished for i=$i. Running Analysis..."
  # Pipe B_BIN output to a file
  "$B_BIN" \
    -bg $base_graph \
    -sg "./../data/daniel-25.fmi" \
    -csv "../../data/histogram_$i.csv" \
    -p 0.0 -k $i | tee -a "./../data/analysis_output.txt"

  echo "Analysis output saved to ../../data/analysis_output.txt"
done

echo "All done."