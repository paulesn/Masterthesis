#!/usr/bin/env bash
set -euo pipefail

# Change these if your binaries live elsewhere
A_BIN="${A_BIN:-./src/daniel/theta-graph/cmake-build-release/createThetaGraph}"   # C++ program A
B_BIN="${B_BIN:-./cmake-build-release/analyse}"   # C++ program B

# The sequence of i values
nums=(24 50 64 100 124)

"ls" "../data"

for i in "${nums[@]}"; do
  echo "Running Daniels Theta spanner generation with theta =$i..."
  "$A_BIN" "-g" "./../data/coastlines-mercator.txt.pruned.wc.txt.graph" "-o" "./../data/daniel.fmi"  "-n" "$i"          # runs A synchronously; waits until it finishes
  echo "A finished for i=$i. Running Analysis..."
  "$B_BIN" "-bg"  "./../data/coastlines-mercator.txt.pruned.wc.txt.graph" "-sg" "./../data/daniel.fmi" "-csv" "../../data/histogram_$i.csv"            # runs B after A completes
done

echo "All done."

