#!/bin/bash
# Helper script to import GLMM.R from common paths

set -e

# Define common source paths
COMMON_PATHS=(
  "D:/USDA/TRACE_DM_DEC/STAT/GLMM.R"
  "$HOME/USDA/TRACE_DM_DEC/STAT/GLMM.R"
  "/mnt/d/USDA/TRACE_DM_DEC/STAT/GLMM.R"
  "./GLMM.R"
  "../GLMM.R"
)

# Target destination
TARGET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_FILE="$TARGET_DIR/GLMM.R"

echo "Searching for GLMM.R in common locations..."
echo "Target: $TARGET_FILE"
echo ""

FOUND=0
for path in "${COMMON_PATHS[@]}"; do
  if [ -f "$path" ]; then
    echo "✓ Found GLMM.R at: $path"
    cp "$path" "$TARGET_FILE"
    echo "✓ Copied to: $TARGET_FILE"
    echo "✓ File size: $(wc -c < "$TARGET_FILE") bytes"
    FOUND=1
    break
  fi
done

if [ $FOUND -eq 0 ]; then
  echo "✗ GLMM.R not found in any common location."
  echo ""
  echo "Checked paths:"
  for path in "${COMMON_PATHS[@]}"; do
    echo "  - $path"
  done
  exit 1
fi

echo ""
echo "Import complete!"
