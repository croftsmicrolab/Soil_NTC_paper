# GNAT domain coordinates for 18 sequenced candidates plus SatA and SatB

The 18 experimentally sequenced protein sequences and the SatA/SatB tree-reference proteins were scanned together against the local Pfam-A database with HMMER 3.4 `hmmscan --cut_ga`.
Overlapping Pfam GNAT-family models were consolidated into physical domain regions. Nineteen proteins contain one predicted GNAT region; HIT_000035 contains two non-overlapping predicted GNAT regions.

Coordinate definitions:

- `aa_alignment_start/end`: residues aligned to the Pfam HMM; use these as the primary GNAT-domain coordinates.
- `aa_envelope_start/end`: the wider HMMER domain envelope; useful as a conservative boundary.
- `cds_nt_*`: 1-based coordinates within the selected true candidate CDS.
- `sequenced_insert_nt_*`: 1-based coordinates within the original sequenced candidate insert, accounting for user-selected offsets such as nt47 or nt9.
- SatA and SatB were supplied as protein-only references, so their nucleotide-coordinate fields are `NA` rather than inferred.
- Coordinates are profile-based predictions, not experimentally determined structural boundaries.

| Sequence | Protein length (aa) | Best Pfam GNAT model | Protein alignment (aa) | Protein envelope (aa) | Selected CDS alignment (nt) | Sequenced insert alignment (nt) | i-Evalue |
|---|---:|---|---:|---:|---:|---:|---:|
| HIT_000122 | 138 | Acetyltransf_1 (PF00583.32) | 23–127 | 6–127 | 67–381 | 67–381 | 5.6e-19 |
| HIT_000135 | 137 | Acetyltransf_1 (PF00583.32) | 4–109 | 1–109 | 10–327 | 10–327 | 3.1e-18 |
| HIT_000242 | 177 | Acetyltransf_1 (PF00583.32) | 50–149 | 39–149 | 148–447 | 148–447 | 2.2e-23 |
| HIT_000091 | 191 | Acetyltransf_1 (PF00583.32) | 49–149 | 40–149 | 145–447 | 145–447 | 8e-22 |
| HIT_000259 | 178 | Acetyltransf_1 (PF00583.32) | 56–150 | 23–150 | 166–450 | 199–483 | 5.6e-23 |
| HIT_000197 | 180 | Acetyltransf_1 (PF00583.32) | 49–149 | 39–149 | 145–447 | 145–447 | 4.6e-22 |
| HIT_000008 | 187 | Acetyltransf_1 (PF00583.32) | 63–159 | 41–159 | 187–477 | 187–477 | 2.3e-20 |
| HIT_000057 | 185 | Acetyltransf_1 (PF00583.32) | 48–155 | 41–155 | 142–465 | 142–465 | 6.6e-21 |
| HIT_000118 | 181 | Acetyltransf_1 (PF00583.32) | 55–153 | 38–153 | 163–459 | 163–459 | 1.3e-19 |
| HIT_000040 | 178 | Acetyltransf_1 (PF00583.32) | 49–148 | 30–148 | 145–444 | 145–444 | 1.8e-22 |
| HIT_000156 | 177 | Acetyltransf_1 (PF00583.32) | 48–149 | 37–149 | 142–447 | 142–447 | 8.9e-24 |
| HIT_000262 | 178 | Acetyltransf_1 (PF00583.32) | 45–151 | 33–151 | 133–453 | 179–499 | 4.3e-21 |
| HIT_000178 | 323 | Acetyltransf_1 (PF00583.32) | 190–310 | 183–310 | 568–930 | 568–930 | 5.3e-10 |
| HIT_000185 | 199 | Acetyltransf_1 (PF00583.32) | 109–197 | 94–199 | 325–591 | 371–637 | 1.9e-11 |
| HIT_000151 | 188 | Acetyltransf_3 (PF13302.14) | 11–152 | 11–152 | 31–456 | 31–456 | 3.2e-34 |
| HIT_000035 region 1 | 307 | Acetyltransf_1 (PF00583.32) | 56–133 | 28–133 | 166–399 | 166–399 | 0.00013 |
| HIT_000035 region 2 | 307 | Acetyltransf_1 (PF00583.32) | 199–296 | 173–296 | 595–888 | 595–888 | 1.2e-13 |
| HIT_000196 | 174 | Acetyltransf_1 (PF00583.32) | 49–148 | 24–149 | 145–444 | 153–452 | 9.8e-16 |
| HIT_000213 | 159 | Acetyltransf_1 (PF00583.32) | 26–130 | 17–130 | 76–390 | 76–390 | 3.6e-24 |
| SatA | 184 | Acetyltransf_1 (PF00583.32) | 50–154 | 40–154 | NA | NA | 9.7e-25 |
| SatB | 196 | Acetyltransf_1 (PF00583.32) | 72–169 | 56–169 | NA | NA | 3.4e-23 |

Files:

- `gnat_domain_coordinates_20.tsv`: all distinct GNAT regions (21 regions across 20 proteins).
- `gnat_primary_domain_coordinates_20.tsv`: one strongest GNAT-domain call per protein.
- `gnat_all_pfam_domain_hits_20.tsv`: every overlapping GA-threshold GNAT-family Pfam hit.
- `pfam_selected_true_18_plus_SatA_SatB.domtblout`: raw HMMER domain table.
- `selected_true_18_plus_SatA_SatB.faa`: the exact 20-protein analysis input.
- `reference_SatA_SatB.faa`: SatA and SatB extracted from the 88-reference tree FASTA.
- `sequence_metadata_20.tsv`: candidate/reference source metadata.
- `run_gnat_domain_mapping.R`: reproducible analysis script.
