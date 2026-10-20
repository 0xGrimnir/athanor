#!/bin/bash
# Report tool versions and command-line interfaces from inside the pack image.
B=/opt/athanor/bio/bin
echo "== versions"
$B/bwa 2>&1 | grep -i version
$B/samtools --version | head -1
$B/AdapterRemoval --version 2>&1 | head -1
$B/mapDamage --version 2>&1 | head -1
ls $B | grep -iE "garg|deam|frag|adptsim|art_" | tr '\n' ' '; echo
echo "== gargammel help"
$B/gargammel --help 2>&1 | head -120 || $B/gargammel.pl --help 2>&1 | head -120
echo "== refs"
cat /opt/athanor/refs/manifest.json | head -60
