# Negative genomic adapter controls

Two 10,000-base windows downloaded on 2026-10-15 from NCBI Nucleotide. FASTA
headers retain source accession and coordinates. Coordinates are one-based,
inclusive; no bases were changed. Tests use these checked-in files without
network access, testing 150-base windows on both strands against synthetic
adapters. These controls supplement randomized and low-complexity negative reads;
they do not establish false-positive rates for a particular sequencing library.

- [E. coli K-12 MG1655 U00096.3, 1–10,000](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U00096.3&seq_start=1&seq_stop=10000&rettype=fasta&retmode=text)
- [Human GRCh38.p14 chromosome 22 NC_000022.11, 20,000,001–20,010,000](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000022.11&seq_start=20000001&seq_stop=20010000&rettype=fasta&retmode=text)
