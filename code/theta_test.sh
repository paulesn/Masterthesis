#!/usr/bin/env bash
set -euo pipefail

# Change these if your binaries live elsewhere
A_BIN="${A_BIN:-../theta-graph/cmake-build-release/createThetaGraph}"   # C++ program A
B_BIN="${B_BIN:-./cmake-build-release/analyse}"   # C++ program B

# The sequence of i values
nums=(12 24 30 40 50 64 75 90 100 111 124)

"ls" "../data"

for i in "${nums[@]}"; do
  echo "Running Daniels Theta spanner generation with theta =$i..."
  "$A_BIN" "-g" "./../data/coastlines-mercator.txt.pruned.wc.txt.shrunk.0025POI.32.txt.graph" "-o" "./../data/daniel.fmi"  "-n" "$i"          # runs A synchronously; waits until it finishes
  echo "A finished for i=$i. Running Analysis..."
  "$B_BIN" "-bg" "../data/0025.32.fmi" "-sg" "./../data/daniel.fmi" "-csv" "../../data/histogram_$i.csv"            # runs B after A completes
done

echo "All done."