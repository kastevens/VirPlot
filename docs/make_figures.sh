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
virplot -g examples/grbv.gff3 -d examples/grbv.sam -y examples/grbv_ictv.yml -f svg --title --smooth --layout circular -o "$out" --name grbv_circular_ictv
virplot -g examples/pvy.gff3  -d examples/pvy.sam  -y examples/pvy.yml  -f svg --title --smooth --grid -o "$out" --name pvy_linear
# Segmented genome: one figure per segment (largest first, bare x-axis except
# the bottom one), then composed into the ICTV stacked figure. The per-segment
# files are scratch; only the stacked figure is kept.
virplot -g examples/tswv.gff3 -d examples/tswv.sam -y examples/tswv.yml -f svg --title --smooth --bare-x -o "$out" --name tswv_segment
python3 bin/stack_figures.py -i "$out"/tswv_segment.*.svg -o "$out"/tswv_stacked.svg
rm -f "$out"/tswv_segment.*.svg
