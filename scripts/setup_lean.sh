#!/usr/bin/env bash
# Script to install Lean 4 version manager (elan) and build the formalization
set -e

echo "=================================================================="
echo "  Lean 4 Formalization Setup for Navier-Stokes & Euler Breakdown"
echo "=================================================================="

# Check if elan is already installed
if command -v elan &> /dev/null; then
    echo "✓ elan is already installed: $(elan --version)"
else
    echo "Installing elan (Lean version manager) via official installer..."
    curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y
    source "$HOME/.elan/env"
    echo "✓ elan successfully installed: $(elan --version)"
fi

# Ensure ~/.elan/bin is in PATH
export PATH="$HOME/.elan/bin:$PATH"

echo ""
echo "Checking Lean toolchain specified in lean-toolchain..."
elan which lean
elan which lake

echo ""
echo "Fetching Mathlib pre-compiled cache (saves hours of compilation time)..."
lake exe cache get || echo "Note: Mathlib cache fetch had warnings; proceeding with lake build..."

echo ""
echo "Building formalization targets..."
lake build NavierStokes
lake build Euler

echo ""
echo "=================================================================="
echo "✓ Formalization build complete!"
echo "To inspect proofs:"
echo "  lake build NavierStokes.ComparatorSolution"
echo "=================================================================="
