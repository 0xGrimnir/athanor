#!/bin/bash
# Report tool versions and command-line interfaces from inside the pack image.
export PATH=/opt/athanor/bio/bin:$PATH
B=/opt/athanor/bio/bin
echo "== binaries"; ls $B | grep -iE "adapter|garg|deam|frag|adpt|art_illumina|mapdamage|minimap|bwa|samtools|perl$" | tr '\n' ' '; echo
echo "== versions"
bwa 2>&1 | grep -i version; samtools --version | head -1; minimap2 --version
for a in AdapterRemoval adapterremoval adapterremoval3; do command -v $a >/dev/null && { echo "$a:"; $a --version 2>&1 | head -2; }; done
mapDamage --version 2>&1 | head -1
echo "== gargammel --help"; gargammel --help 2>&1 | head -150
echo "== fragSim"; fragSim 2>&1 | head -60
echo "== deamSim"; deamSim 2>&1 | head -60
echo "== AdapterRemoval --help"; (AdapterRemoval --help || adapterremoval3 --help) 2>&1 | head -80
echo "== mapDamage --help"; mapDamage --help 2>&1 | head -120
echo "== licenses (conda package metadata)"
for f in /opt/athanor/bio/conda-meta/*.json; do python3 -c "import json,sys;d=json.load(open(sys.argv[1]));print(d['name'],d['version'],'|',d.get('license','?'))" "$f"; done | grep -iE "^(gargammel|art|bwa|samtools|htslib|adapterremoval|mapdamage2|minimap2|r-base|perl|python) " 
