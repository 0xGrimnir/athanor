Athanor 0.3.0: the Ancient DNA Lab.

- **New tab: Ancient DNA.** Simulate ancient reads from a real mitochondrial genome (woolly mammoth by default), map them to relatives' genomes (African and Asian elephants), and watch the reconstruction: damage patterns recovered by mapDamage2, a consensus where every base carries a probability, checked against the true sequence, and iterative mapping that removes most of the reference bias.
- Uses the field's own tools: gargammel, AdapterRemoval, BWA, mapDamage2, minimap2.
- **Requires Linux engines pack 2** (about 2.4 GB). Click **Install Linux engines** in the top bar after updating; it replaces the previous pack. No GPU is needed for the ancient-DNA tools.
