#!/usr/bin/env bash

# build.sh ────────────────────────────────────────────────────────────────
# Simple Linux build script for the atom visualizer (raylib)

set -euo pipefail

# Configuration ────────────────────────────────────────────────────────────

CC="gcc"
CFLAGS="-O2 -s -Wall -Wextra"
LIBS="-lraylib -lGL -lm -lpthread -ldl -lX11"
# LIBS="-lraylib"   # ← try this first if raylib is installed system-wide

OUTPUT="atom"
SOURCE="atom.c"

# Optional: if you have local raylib (uncomment and adjust path)
# LOCAL_RAYLIB="-I./raylib/src ./raylib/src/libraylib.a"

# ──────────────────────────────────────────────────────────────────────────

echo "Building ${OUTPUT} ..."

${CC} ${SOURCE} -o ${OUTPUT} \
    ${CFLAGS} \
    ${LOCAL_RAYLIB:-} \
    ${LIBS}
sudo apt update && sudo apt install upx-ucl
upx --ultra-brute --lzma ${OUTPUT}
if [ $? -eq 0 ]; then
    echo ""
    echo "Build successful!"
    echo ""
    echo "Run examples:"
    echo "  ./atom 1 30     # Hydrogen"
    echo "  ./atom 8 60     # Oxygen"
    echo "  ./atom 26 40    # Iron"
    echo "  ./atom 118 80   # Oganesson (very large zoom recommended)"
    echo "runnig example with 26 protons and 40x zoom..."
    # Uncomment if you want to run automatically:
    ./atom 26 40
else
    echo ""
    echo "Build failed."
    echo "Common fixes:"
    echo "  • Install raylib:  sudo apt install libraylib-dev"
    echo "  • Missing libs:    sudo apt install libgl1-mesa-dev libx11-dev ..."
    echo "  • Try simpler link:  ${CC} ${SOURCE} -o ${OUTPUT} -lraylib"
fi