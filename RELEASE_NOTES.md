Athanor 0.9.0: DNA Structure, complete.

- The helix's backbone is now closed with a real force field (OpenMM, AMBER OL15). The helix can also move: molecular dynamics at the temperature and salt you choose, on your graphics card when it can.
- The shape DNA takes in solution, predicted from its sequence (Deep DNAshape), on a plot and as colour on the helix.
- What holds the helix together: the pairing energy of every base pair and the stacking energy of every step, with a DFT check for any pair you pick.
- Ancient damage on a molecule: single-stranded ends and the cytosines turned to uracil, from the Briggs model or from the damage your Ancient DNA Lab run measured.
- A zoom through the scales:
  - your sequence wrapped on a nucleosome;
  - proteins bending DNA (TBP, CAP, IHF and others);
  - a six-nucleosome chromatin fibre;
  - a whole human nucleus, where a Personal Genome variant can be found.
- Every new feature is checked against published data or the original programs; the README lists the results.
