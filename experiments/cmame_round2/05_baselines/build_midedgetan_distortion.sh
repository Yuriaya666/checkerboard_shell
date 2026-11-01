#!/usr/bin/env bash
set -euo pipefail

source_root="${LIBSHELL_SOURCE_ROOT:-/tmp/platepinn_checkerboard_benchmark_sources/libshell}"
eigen_root="${EIGEN_SOURCE_ROOT:-/tmp/platepinn_checkerboard_benchmark_sources/eigen}"
output="${1:-/tmp/checkerboard_round2_midedgetan_distortion}"
driver_root="$(cd "$(dirname "$0")" && pwd)"

g++ -std=c++17 -O3 -DNDEBUG \
  -I"$source_root/include" -I"$source_root/src" -I"$eigen_root" \
  "$driver_root/midedgetan_distortion_driver.cpp" \
  "$source_root/src/ElasticShell.cpp" \
  "$source_root/src/GeometryDerivatives.cpp" \
  "$source_root/src/MeshConnectivity.cpp" \
  "$source_root/src/MaterialModel/StVKMaterial.cpp" \
  "$source_root/src/SecondFundamentalForm/MidedgeAverageFormulation.cpp" \
  "$source_root/src/SecondFundamentalForm/MidedgeAngleTanFormulation.cpp" \
  "$source_root/src/SecondFundamentalForm/MidedgeAngleSinFormulation.cpp" \
  "$source_root/src/SecondFundamentalForm/MidedgeAngleThetaFormulation.cpp" \
  -o "$output"

echo "$output"
