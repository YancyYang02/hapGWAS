#!/bin/bash
# End-to-end example on the toy dataset (see README, "Example").
set -e
cd "$(dirname "$0")"
python make_toy.py
hapdc prep toy/panel.vcf.gz auto toy/height.txt run/prep
hapdc ms   run/prep run/ms 0 --covar toy/covar.txt
echo "planted signal:"; cat toy/truth.txt
echo "significant anchors (Bonferroni):"
awk -F'\t' 'NR==1 || $10==1' run/ms/height.anchors.tsv
