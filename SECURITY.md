# Security and privacy

## Reporting a vulnerability

Please report security issues privately through GitHub's **Report a vulnerability** button on the repository's Security tab, rather than in a public issue. Include steps to reproduce and the version affected.

## Privacy design

Athanor handles some sensitive inputs, personal genetic reports above all, so it follows these rules:

- **Local processing.** DNA reports and genotype files are parsed on your machine. They are never uploaded, sent to a server, or included in telemetry.
- **No telemetry.** The app sends no usage data.
- **Network use is explicit.** The app contacts the network only to check for updates and to fetch public reference data you request (for example Ensembl, NCBI, PDB, AlphaFold, ClinVar or gnomAD). Your own files are never part of those requests.
- **Signed updates.** Updates are cryptographically signed, and the app refuses any update whose signature doesn't verify.
- **The repository excludes genetic data.** `.gitignore` blocks common consumer-DNA and sequencing file patterns so personal data can't be committed by accident.
