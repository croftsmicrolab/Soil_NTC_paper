#!/bin/sh
set -eu
package_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

printf "%s\n" "Running Figures/AF3_structure/RMSD.R"
(cd "$package_dir/Figures/AF3_structure" && Rscript --vanilla RMSD.R)

printf "%s\n" "Running Figures/AF3_structure/Reproducibility/annotation_source/run_gnat_domain_mapping.R"
(cd "$package_dir/Figures/AF3_structure/Reproducibility/annotation_source" && Rscript --vanilla run_gnat_domain_mapping.R)

printf "%s\n" "Running Figures/Focused_tree/Focused_tree.R"
(cd "$package_dir/Figures/Focused_tree" && Rscript --vanilla Focused_tree.R)

printf "%s\n" "Running Figures/GO/plot_gene_ontology_analysis.R"
(cd "$package_dir/Figures/GO" && Rscript --vanilla plot_gene_ontology_analysis.R)

printf "%s\n" "Running Figures/Microdilution/KAN.R"
(cd "$package_dir/Figures/Microdilution" && Rscript --vanilla KAN.R)

printf "%s\n" "Running Figures/Microdilution/NTC.R"
(cd "$package_dir/Figures/Microdilution" && Rscript --vanilla NTC.R)

printf "%s\n" "Running Figures/Microdilution/SPEC.R"
(cd "$package_dir/Figures/Microdilution" && Rscript --vanilla SPEC.R)

printf "%s\n" "Running Figures/Microdilution/ami.R"
(cd "$package_dir/Figures/Microdilution" && Rscript --vanilla ami.R)

printf "%s\n" "Running Figures/Microdilution/gen.R"
(cd "$package_dir/Figures/Microdilution" && Rscript --vanilla gen.R)

printf "%s\n" "Running Figures/Overview_tree/overview_tree.R"
(cd "$package_dir/Figures/Overview_tree" && Rscript --vanilla overview_tree.R)

printf "%s\n" "Running Figures/Read_length_distribution/plot_read_length_distribution.R"
(cd "$package_dir/Figures/Read_length_distribution" && Rscript --vanilla plot_read_length_distribution.R)

printf "%s\n" "Running Figures/Similarity_heatmap/similarity_heatmap.R"
(cd "$package_dir/Figures/Similarity_heatmap" && Rscript --vanilla similarity_heatmap.R)

printf "%s\n" "Running Figures/orf_cluster/plot_clone24_reporf_clusters.R"
(cd "$package_dir/Figures/orf_cluster" && Rscript --vanilla plot_clone24_reporf_clusters.R)

printf "%s\n" "Running Figures/orf_cluster/plot_orf_cluster_curve.R"
(cd "$package_dir/Figures/orf_cluster" && Rscript --vanilla plot_orf_cluster_curve.R)

printf "%s\n" "Running Figures/raw_sequence_cluster/plot_clone24_raw_sequence_clusters.R"
(cd "$package_dir/Figures/raw_sequence_cluster" && Rscript --vanilla plot_clone24_raw_sequence_clusters.R)

printf "%s\n" "Running Figures/raw_sequence_cluster/plot_raw_sequence_cluster_curve.R"
(cd "$package_dir/Figures/raw_sequence_cluster" && Rscript --vanilla plot_raw_sequence_cluster_curve.R)

printf "%s\n" "Running Figures/taxonomy/taxonomy_plot.R"
(cd "$package_dir/Figures/taxonomy" && Rscript --vanilla taxonomy_plot.R)

printf "%s\n" "Running Figures/test_strip_heatmap/test_strip_log2fc_heatmap.R"
(cd "$package_dir/Figures/test_strip_heatmap" && Rscript --vanilla test_strip_log2fc_heatmap.R)
