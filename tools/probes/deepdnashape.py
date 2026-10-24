"""Probe: Deep DNAshape checkpoints (variable names/shapes) and a reference prediction."""
import os, subprocess, sys
subprocess.run(["git", "clone", "-q", "--depth", "1", "https://github.com/JinsenLi/deepDNAshape.git", "/tmp/dds"], check=True)
print(subprocess.run(["git", "-C", "/tmp/dds", "log", "-1", "--format=%H %ci"], capture_output=True, text=True).stdout)
import tensorflow as tf
print("tf", tf.__version__)
for f in ("MGW", "HelT", "MGW-FL"):
    r = tf.train.load_checkpoint(f"/tmp/dds/deepDNAshape/models/{f}")
    for k, s in sorted(r.get_variable_to_shape_map().items()):
        if "OPTIMIZER" in k or "optimizer" in k: continue
        print(f, k, s)
sys.path.insert(0, "/tmp/dds")
from deepDNAshape import predictor
p = predictor.predictor()
for f in ("MGW", "HelT", "Roll", "ProT"):
    print(f, p.predict(f, "CGCGAATTCGCG", 4).round(4).tolist())
