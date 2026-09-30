#!/usr/bin/env sh
# Regenerate the figures embedded in the documentation (README, docs/).
# Run from the repository root:  sh docs/make_figures.sh
set -e
out=docs/figures
virplot -g examples/byv.gff3  -d examples/byv.sam  -y examples/byv.yml  -f svg --title -o "$out" --name byv_linear
virplot -g examples/grbv.gff3 -d examples/grbv.sam -y examples/grbv.yml -f svg --title --legend -o "$out" --name grbv_linear
virplot -g examples/grbv.gff3 -d examples/grbv.sam -y examples/grbv.yml -f svg --title --legend --layout circular -o "$out" --name grbv_circular
virplot -g examples/sample.gff3 -d examples/sample.dep -y examples/spec.yml -f svg --title --legend --smooth --shade-breaks --grid -o "$out" --name sample_linear
