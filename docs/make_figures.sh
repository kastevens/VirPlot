#!/usr/bin/env sh
# Regenerate the figures embedded in the documentation (README, docs/).
# Run from the repository root:  sh docs/make_figures.sh
set -e
out=docs/figures
rm -f "$out"/*.svg          # virplot never overwrites, so clear the old set first
virplot -g examples/byv.gff3  -d examples/byv.sam  -y examples/byv.yml  -f svg --title --smooth --grid -o "$out" --name byv_linear
virplot -g examples/grbv.gff3 -d examples/grbv.sam -y examples/grbv.yml -f svg --title --legend --smooth --grid -o "$out" --name grbv_linear
virplot -g examples/grbv.gff3 -d examples/grbv.sam -y examples/grbv.yml -f svg --title --legend --smooth --layout circular -o "$out" --name grbv_circular
virplot -g examples/sample.gff3 -d examples/sample.dep -y examples/spec.yml -f svg --title --legend --smooth --shade-breaks --grid -o "$out" --name sample_linear
