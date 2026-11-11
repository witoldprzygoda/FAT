#!/bin/bash
# rebuild.sh -- builds the deck, applies the theme, validates it and renders every slide.
#   ./rebuild.sh          deck only
#   ./rebuild.sh assets   also figures, equations and table numbers (needs data/deck_numbers.json)
# Environment: TOOLS = directory with node_modules (pptxgenjs 4.0.1); SKILL = directory of the pptx skill;
#              RENDER = output directory of the slide images (default: ./render).
set -e
cd "$(dirname "$(readlink -f "$0")")"
: "${NODE_BIN:=/cvmfs/sft.cern.ch/lcg/releases/nodejs/22.14.0-07111/x86_64-el9-gcc13-opt/bin}"
: "${TOOLS:?set TOOLS to the directory with node_modules}"
: "${SKILL:?set SKILL to the pptx skill directory}"
: "${RENDER:=render}"
export PATH="$NODE_BIN:$PATH" NODE_PATH="$TOOLS/node_modules" SKILL
DECK=Tlo_kombinatoryczne_SMASH_PL
if [ "$1" = assets ]; then python3 make_assets.py >/dev/null; fi
node build_deck.js
PYTHONPATH= python3 "$SKILL/scripts/office/validate.py" "$DECK.pptx" || true
mkdir -p "$RENDER"; rm -f "$RENDER"/s-*.png "$RENDER/$DECK.pdf"
python3 "$SKILL/scripts/office/soffice.py" --headless --convert-to pdf --outdir "$RENDER" "$DECK.pptx" >/dev/null 2>&1
mutool draw -r 80 -o "$RENDER/s-%02d.png" "$RENDER/$DECK.pdf" >/dev/null 2>&1
ls "$RENDER" | wc -l
