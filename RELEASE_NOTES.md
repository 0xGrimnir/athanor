Athanor 0.10.0: Deep Time begins — ancient nuclear genomes.

- A new tab reads an ancient genome as it is deposited: a BAM aligned to GRCh37, on your own computer, nothing uploaded. A whole-genome file can be hundreds of gigabytes; only the positions needed are read from it.
- What the alignment is: which reference build, how many reads on each chromosome (from the index, so it is quick), the libraries and programs the file records, and how long the molecules are.
- Is this DNA ancient? Deamination towards both ends of the molecules, with cytosines followed by guanine counted separately, because methylated CpG deaminates several times faster. Give the reference genome the reads were mapped to and every aligned base is compared with it, as mapDamage2 does.
- Genetic sex from the depth at the reference panel's own positions, which reads capture data and shotgun data alike.
- Genotypes at the panel's sites — pseudo-haploid, as published ancient panels are built, or diploid where coverage allows — and from them the same reports the Personal Genome tab makes: ancient ancestry, principal components and closest ancient individuals, Neanderthal and Denisovan ancestry, maternal and paternal lines.
- Heterozygosity over a region, counted from the reads, with a study's own callable-regions mask; the transversion rate is shown beside the total, and the app says when deamination is inflating the total.
- Checked against published ancient genomes: calls agree with the Allen Ancient DNA Resource's own genotypes for four 1240K individuals at 98.3–99.9% of shared sites, sex matches all four, and in two megabases of Vindija 33.19 the diploid calls agree with that study's call set at 99.86% of 1.3 million sites.
- Ancient DNA Lab: the simulator is now gargammel 1.1.5, built from source, whose damage model is the published one. All three Briggs parameters now come back from the simulated reads; δs used to read about twice its true value.
- This release needs Linux engines pack 3; the app asks to install it on first run.
