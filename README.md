# Soil NTC paper

Code, input data and final illustrations for the Soil NTC manuscript.
Each analysis folder contains the files needed by its scripts.

## Figure guide

| Figure | Folder under `Figures/` |
|---|---|
| 1 | `SAT_figure` |
| 2 | `functional_metagenomic_pipline_TSC_version` |
| 3 | `pipeline` |
| 4 | `Overview_tree` |
| 5 | `Gene schema`, `taxonomy` |
| 6 | `Agar_dilution`, `Focused_tree` |
| 7 | `Microdilution` |
| 8 | `AF3_structure` |
| S1 | `Selection_plate`, `Read_length_distribution` |
| S2 | `raw_sequence_cluster` |
| S3 | `orf_cluster` |
| S4 | `taxonomy` |
| S5 | `GO` |
| S6 | `Similarity_heatmap` |
| S7 | `Microdilution` |
| S8 | `sequence_analysis` |

`STC_structure`, `test_strip`, `test_strip_heatmap` and the gentamicin files in
`Microdilution` contain additional original assets or analyses.
Folders without plotting code contain final figure exports. Scripts reproduce
individual panels; composite publication layouts were assembled manually.

## Requirements

Tested with R 4.6.0 and UCSF ChimeraX 1.12 on macOS.
R package names and versions are listed in `environment/R_package_versions.tsv`.
`Biostrings`, `ComplexHeatmap` and `ggtree` are Bioconductor packages; the remaining
listed packages are from CRAN or included with R.

Structure-image composition requires Python with NumPy and Pillow:

```sh
python3 -m pip install -r environment/requirements.txt
```

The structure-composition script uses the standard macOS Arial font paths.
On other systems, set its two font paths to the installed Arial font files.

## Run

Each command block below starts from the repository root.

Run all R analyses from the repository root:

```sh
sh run_R_figures.sh
```

To run an individual script, change into its folder first, for example:

```sh
cd Figures/Overview_tree
Rscript --vanilla overview_tree.R
```

Generate structural alignments and images:

```sh
cd Figures/AF3_structure/Reproducibility
sh run_structure.sh
```

If needed, set `CHIMERAX_BIN` to the ChimeraX executable and `FIGURE_PYTHON` to
Python with NumPy/Pillow. Rendering requires graphics support. `RMSD.R` in the
parent folder draws the RMSD plot from the included `RMSD.txt`.

The GNAT annotation script uses the included original Pfam scan by default.
An optional rescan requires HMMER 3.4 and Pfam 38.2:

```sh
cd Figures/AF3_structure/Reproducibility/annotation_source
PFAM_HMM=/path/to/Pfam-A.hmm Rscript --vanilla run_gnat_domain_mapping.R --rescan
```

Figures and generated tables are written beside the scripts or into their
output subdirectories. Reproduction starts from the included processed inputs;
upstream sequencing analysis and AlphaFold prediction are not rerun.
