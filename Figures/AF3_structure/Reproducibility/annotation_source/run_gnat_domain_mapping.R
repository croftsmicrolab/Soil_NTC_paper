#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Run from this annotation_source directory. Reuse the archived Pfam scan by
# default; --rescan requires HMMER on PATH and PFAM_HMM set to Pfam-A.hmm.
input_dir <- normalizePath(".", mustWork = TRUE)
output_dir <- file.path(input_dir, "generated")
input_faa <- file.path(input_dir, "selected_true_proteins_18.faa")
input_fna <- file.path(input_dir, "selected_true_regions_18.fna")
input_qc <- file.path(input_dir, "selected_true_regions_18_qc.tsv")
input_reference_faa <- file.path(input_dir, "crofts_and_paper_method3_references_88.faa")
input_reference_mapping <- file.path(input_dir, "reference_label_mapping_88.tsv")
rescan <- "--rescan" %in% commandArgs(trailingOnly = TRUE)
hmmscan_bin <- Sys.which("hmmscan")
pfam_hmm <- Sys.getenv("PFAM_HMM", unset = "")
archived_domtblout <- file.path(input_dir, "pfam_selected_true_18_plus_SatA_SatB.domtblout")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required_files <- c(input_faa, input_fna, input_qc, input_reference_faa, input_reference_mapping,
                    if (rescan) c(hmmscan_bin, pfam_hmm) else archived_domtblout)
missing_files <- required_files[!nzchar(required_files) | !file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files (for --rescan, set PFAM_HMM and put hmmscan on PATH): ",
       paste(missing_files, collapse = ", "))
}

protein_copy <- file.path(output_dir, "selected_true_proteins_18.faa")
nucleotide_copy <- file.path(output_dir, "selected_true_regions_18.fna")
qc_copy <- file.path(output_dir, "selected_true_regions_18_qc.tsv")
reference_copy <- file.path(output_dir, "reference_SatA_SatB.faa")
combined_faa <- file.path(output_dir, "selected_true_18_plus_SatA_SatB.faa")
invisible(file.copy(input_faa, protein_copy, overwrite = TRUE))
invisible(file.copy(input_fna, nucleotide_copy, overwrite = TRUE))
invisible(file.copy(input_qc, qc_copy, overwrite = TRUE))

read_fasta <- function(path) {
  lines <- readLines(path)
  header_positions <- which(grepl("^>", lines))
  if (length(header_positions) == 0L) stop("No FASTA records in ", path)
  end_positions <- c(header_positions[-1] - 1L, length(lines))
  sequences <- vapply(seq_along(header_positions), function(i) {
    paste(lines[(header_positions[i] + 1L):end_positions[i]], collapse = "")
  }, character(1))
  names(sequences) <- sub("^>", "", lines[header_positions])
  sequences
}

write_fasta <- function(sequences, path, width = 60L) {
  output <- unlist(lapply(names(sequences), function(id) {
    sequence <- sequences[[id]]
    starts <- seq.int(1L, nchar(sequence), by = width)
    chunks <- substring(sequence, starts, pmin(starts + width - 1L, nchar(sequence)))
    c(paste0(">", id), chunks)
  }), use.names = FALSE)
  writeLines(output, path)
}

candidate_proteins <- read_fasta(protein_copy)
all_references <- read_fasta(input_reference_faa)
reference_ids <- c("SatA", "SatB")
if (!all(reference_ids %in% names(all_references))) {
  stop("SatA and/or SatB are absent from the tree reference FASTA")
}
selected_references <- all_references[reference_ids]
write_fasta(selected_references, reference_copy)
combined_proteins <- c(candidate_proteins, selected_references)
if (length(combined_proteins) != 20L || anyDuplicated(names(combined_proteins))) {
  stop("Expected 20 unique proteins after adding SatA and SatB")
}
write_fasta(combined_proteins, combined_faa)

fasta_headers <- sub("^>", "", grep("^>", readLines(combined_faa), value = TRUE))
if (length(fasta_headers) != 20L || anyDuplicated(fasta_headers)) {
  stop("Expected 20 unique protein FASTA records; found ", length(fasta_headers))
}

domtblout <- file.path(output_dir, "pfam_selected_true_18_plus_SatA_SatB.domtblout")
hmmscan_text <- file.path(output_dir, "pfam_selected_true_18_plus_SatA_SatB.hmmscan.txt")
hmmscan_stderr <- file.path(output_dir, "pfam_selected_true_18_plus_SatA_SatB.hmmscan.stderr.log")

hmmscan_args <- c(
  "--cpu", "8",
  "--cut_ga",
  "--domtblout", domtblout,
  pfam_hmm,
  combined_faa
)

if (rescan) {
  status <- system2(hmmscan_bin, args = hmmscan_args,
                    stdout = hmmscan_text, stderr = hmmscan_stderr)
  if (!identical(status, 0L)) {
    stop("hmmscan failed with exit status ", status, "; see ", hmmscan_stderr)
  }
} else {
  # Confirm that the regenerated scan input is exactly the archived input.
  archived_faa <- file.path(input_dir, "selected_true_18_plus_SatA_SatB.faa")
  if (!identical(unname(tools::md5sum(combined_faa)),
                 unname(tools::md5sum(archived_faa)))) {
    stop("Combined protein FASTA differs from the archived Pfam scan input.")
  }
  if (!file.copy(archived_domtblout, domtblout, overwrite = TRUE)) {
    stop("Could not copy the archived Pfam scan.")
  }
}

parse_domtblout <- function(path) {
  lines <- readLines(path)
  lines <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]
  parsed <- lapply(lines, function(line) {
    fields <- strsplit(trimws(line), "[[:space:]]+")[[1]]
    if (length(fields) < 22L) stop("Malformed domtblout line: ", line)
    c(fields[1:22], paste(fields[-seq_len(22)], collapse = " "))
  })
  mat <- do.call(rbind, parsed)
  column_names <- c(
    "pfam_name", "pfam_accession", "pfam_model_length", "candidate_id", "query_accession",
    "protein_length_aa", "full_evalue", "full_score", "full_bias", "domain_number",
    "domain_count", "domain_c_evalue", "domain_i_evalue", "domain_score", "domain_bias",
    "hmm_start", "hmm_end", "aa_alignment_start", "aa_alignment_end", "aa_envelope_start",
    "aa_envelope_end", "mean_posterior_probability", "pfam_description"
  )
  colnames(mat) <- column_names
  out <- as.data.frame(mat, check.names = FALSE)
  numeric_columns <- c(
    "pfam_model_length", "protein_length_aa", "full_evalue", "full_score", "full_bias",
    "domain_number", "domain_count", "domain_c_evalue", "domain_i_evalue", "domain_score",
    "domain_bias", "hmm_start", "hmm_end", "aa_alignment_start", "aa_alignment_end",
    "aa_envelope_start", "aa_envelope_end", "mean_posterior_probability"
  )
  out[numeric_columns] <- lapply(out[numeric_columns], as.numeric)
  out
}

all_hits <- parse_domtblout(domtblout)
gnat_hits <- all_hits[
  grepl("GNAT", all_hits$pfam_description, ignore.case = TRUE) |
    grepl("^Acetyltransf", all_hits$pfam_name, ignore.case = TRUE),
  , drop = FALSE
]

if (nrow(gnat_hits) == 0L) stop("No GNAT-family Pfam hits were detected")

gnat_hits <- gnat_hits[order(
  gnat_hits$candidate_id,
  gnat_hits$domain_i_evalue,
  -gnat_hits$domain_score,
  gnat_hits$pfam_name
), , drop = FALSE]

# Consolidate overlapping Pfam models into distinct physical GNAT regions.
# Most proteins have several overlapping Acetyltransf_* models describing the
# same fold. HIT_000035 contains two non-overlapping GNAT regions and therefore
# must retain both rather than only the strongest hit.
gnat_hits$gnat_region_number <- NA_integer_
for (id in unique(gnat_hits$candidate_id)) {
  idx <- which(gnat_hits$candidate_id == id)
  ordered_idx <- idx[order(
    gnat_hits$aa_alignment_start[idx],
    gnat_hits$aa_alignment_end[idx],
    gnat_hits$domain_i_evalue[idx]
  )]
  region_number <- 0L
  current_end <- -Inf
  for (j in ordered_idx) {
    if (gnat_hits$aa_alignment_start[j] > current_end) {
      region_number <- region_number + 1L
      current_end <- gnat_hits$aa_alignment_end[j]
    } else {
      current_end <- max(current_end, gnat_hits$aa_alignment_end[j])
    }
    gnat_hits$gnat_region_number[j] <- region_number
  }
}

qc <- read.delim(qc_copy, check.names = FALSE)
required_qc_columns <- c(
  "candidate_id", "construct_label", "round", "sequencing_batch",
  "used_start_1based_inclusive", "selected_cds_nt", "selected_protein_aa"
)
if (!all(required_qc_columns %in% names(qc))) {
  stop("QC table is missing required columns")
}

reference_mapping <- read.delim(input_reference_mapping, check.names = FALSE)
reference_mapping <- reference_mapping[match(reference_ids, reference_mapping$tree_label), , drop = FALSE]
if (anyNA(reference_mapping$tree_label)) stop("SatA/SatB reference metadata were not found")

candidate_metadata <- data.frame(
  candidate_id = qc$candidate_id,
  sequence_type = "Sequenced candidate",
  construct_label = qc$construct_label,
  round = qc$round,
  sequencing_batch = qc$sequencing_batch,
  used_start_1based_inclusive = qc$used_start_1based_inclusive,
  selected_cds_nt = qc$selected_cds_nt,
  selected_protein_aa = qc$selected_protein_aa,
  reference_original_id = NA_character_,
  nucleotide_coordinates_available = TRUE,
  stringsAsFactors = FALSE
)
reference_metadata <- data.frame(
  candidate_id = reference_ids,
  sequence_type = "Reference",
  construct_label = reference_ids,
  round = "Reference",
  sequencing_batch = NA_character_,
  used_start_1based_inclusive = NA_integer_,
  selected_cds_nt = NA_integer_,
  selected_protein_aa = as.integer(nchar(selected_references[reference_ids])),
  reference_original_id = reference_mapping$original_id,
  nucleotide_coordinates_available = FALSE,
  stringsAsFactors = FALSE
)
metadata <- rbind(candidate_metadata, reference_metadata)
write.table(
  metadata,
  file.path(output_dir, "sequence_metadata_20.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

top_index <- !duplicated(gnat_hits$candidate_id)
top_hits <- gnat_hits[top_index, , drop = FALSE]
top_hits <- merge(metadata, top_hits, by = "candidate_id", all.x = TRUE, sort = FALSE)
top_hits <- top_hits[match(metadata$candidate_id, top_hits$candidate_id), , drop = FALSE]

if (anyNA(top_hits$pfam_name)) {
  stop("At least one of the 20 proteins lacks a GA-threshold GNAT-family Pfam hit: ",
       paste(top_hits$candidate_id[is.na(top_hits$pfam_name)], collapse = ", "))
}

top_hits$domain_length_aa <- top_hits$aa_alignment_end - top_hits$aa_alignment_start + 1L
has_nt <- top_hits$nucleotide_coordinates_available
top_hits$cds_nt_alignment_start <- ifelse(has_nt, (top_hits$aa_alignment_start - 1L) * 3L + 1L, NA_integer_)
top_hits$cds_nt_alignment_end <- ifelse(has_nt, top_hits$aa_alignment_end * 3L, NA_integer_)
top_hits$cds_nt_envelope_start <- ifelse(has_nt, (top_hits$aa_envelope_start - 1L) * 3L + 1L, NA_integer_)
top_hits$cds_nt_envelope_end <- ifelse(has_nt, top_hits$aa_envelope_end * 3L, NA_integer_)
top_hits$sequenced_insert_nt_alignment_start <- top_hits$used_start_1based_inclusive + top_hits$cds_nt_alignment_start - 1L
top_hits$sequenced_insert_nt_alignment_end <- top_hits$used_start_1based_inclusive + top_hits$cds_nt_alignment_end - 1L
top_hits$sequenced_insert_nt_envelope_start <- top_hits$used_start_1based_inclusive + top_hits$cds_nt_envelope_start - 1L
top_hits$sequenced_insert_nt_envelope_end <- top_hits$used_start_1based_inclusive + top_hits$cds_nt_envelope_end - 1L

hit_counts <- table(gnat_hits$candidate_id)
top_hits$number_of_overlapping_gnat_pfam_models <- as.integer(hit_counts[top_hits$candidate_id])
top_hits$all_gnat_pfam_models <- vapply(top_hits$candidate_id, function(id) {
  rows <- gnat_hits[gnat_hits$candidate_id == id, , drop = FALSE]
  paste(sprintf("%s (%s)", rows$pfam_name, rows$pfam_accession), collapse = "; ")
}, character(1))

result <- data.frame(
  sequence_id = top_hits$candidate_id,
  sequence_type = top_hits$sequence_type,
  construct_label = top_hits$construct_label,
  round = top_hits$round,
  sequencing_batch = top_hits$sequencing_batch,
  reference_original_id = top_hits$reference_original_id,
  nucleotide_coordinates_available = top_hits$nucleotide_coordinates_available,
  protein_length_aa = top_hits$selected_protein_aa,
  best_gnat_pfam_name = top_hits$pfam_name,
  best_gnat_pfam_accession = top_hits$pfam_accession,
  pfam_description = top_hits$pfam_description,
  aa_alignment_start = top_hits$aa_alignment_start,
  aa_alignment_end = top_hits$aa_alignment_end,
  aa_envelope_start = top_hits$aa_envelope_start,
  aa_envelope_end = top_hits$aa_envelope_end,
  domain_length_aa = top_hits$domain_length_aa,
  cds_nt_alignment_start = top_hits$cds_nt_alignment_start,
  cds_nt_alignment_end = top_hits$cds_nt_alignment_end,
  cds_nt_envelope_start = top_hits$cds_nt_envelope_start,
  cds_nt_envelope_end = top_hits$cds_nt_envelope_end,
  sequenced_insert_nt_alignment_start = top_hits$sequenced_insert_nt_alignment_start,
  sequenced_insert_nt_alignment_end = top_hits$sequenced_insert_nt_alignment_end,
  sequenced_insert_nt_envelope_start = top_hits$sequenced_insert_nt_envelope_start,
  sequenced_insert_nt_envelope_end = top_hits$sequenced_insert_nt_envelope_end,
  domain_i_evalue = top_hits$domain_i_evalue,
  domain_score = top_hits$domain_score,
  mean_posterior_probability = top_hits$mean_posterior_probability,
  number_of_overlapping_gnat_pfam_models = top_hits$number_of_overlapping_gnat_pfam_models,
  all_gnat_pfam_models = top_hits$all_gnat_pfam_models,
  stringsAsFactors = FALSE
)

region_groups <- split(
  seq_len(nrow(gnat_hits)),
  interaction(gnat_hits$candidate_id, gnat_hits$gnat_region_number, drop = TRUE, lex.order = TRUE)
)
region_best_indices <- vapply(region_groups, function(idx) {
  idx[order(gnat_hits$domain_i_evalue[idx], -gnat_hits$domain_score[idx], gnat_hits$pfam_name[idx])][1]
}, integer(1))
region_hits <- gnat_hits[region_best_indices, , drop = FALSE]
region_hits <- merge(metadata, region_hits, by = "candidate_id", all.y = TRUE, sort = FALSE)
region_hits <- region_hits[order(
  match(region_hits$candidate_id, metadata$candidate_id),
  region_hits$gnat_region_number
), , drop = FALSE]

region_hits$gnat_region_count <- as.integer(table(region_hits$candidate_id)[region_hits$candidate_id])
region_hits$is_primary_region <- ave(
  region_hits$domain_i_evalue,
  region_hits$candidate_id,
  FUN = function(x) x == min(x)
) == 1
region_hits$domain_length_aa <- region_hits$aa_alignment_end - region_hits$aa_alignment_start + 1L
has_nt <- region_hits$nucleotide_coordinates_available
region_hits$cds_nt_alignment_start <- ifelse(has_nt, (region_hits$aa_alignment_start - 1L) * 3L + 1L, NA_integer_)
region_hits$cds_nt_alignment_end <- ifelse(has_nt, region_hits$aa_alignment_end * 3L, NA_integer_)
region_hits$cds_nt_envelope_start <- ifelse(has_nt, (region_hits$aa_envelope_start - 1L) * 3L + 1L, NA_integer_)
region_hits$cds_nt_envelope_end <- ifelse(has_nt, region_hits$aa_envelope_end * 3L, NA_integer_)
region_hits$sequenced_insert_nt_alignment_start <- region_hits$used_start_1based_inclusive + region_hits$cds_nt_alignment_start - 1L
region_hits$sequenced_insert_nt_alignment_end <- region_hits$used_start_1based_inclusive + region_hits$cds_nt_alignment_end - 1L
region_hits$sequenced_insert_nt_envelope_start <- region_hits$used_start_1based_inclusive + region_hits$cds_nt_envelope_start - 1L
region_hits$sequenced_insert_nt_envelope_end <- region_hits$used_start_1based_inclusive + region_hits$cds_nt_envelope_end - 1L
region_hits$overlapping_pfam_models <- vapply(seq_len(nrow(region_hits)), function(i) {
  rows <- gnat_hits[
    gnat_hits$candidate_id == region_hits$candidate_id[i] &
      gnat_hits$gnat_region_number == region_hits$gnat_region_number[i],
    , drop = FALSE
  ]
  paste(sprintf("%s (%s)", rows$pfam_name, rows$pfam_accession), collapse = "; ")
}, character(1))

region_result <- data.frame(
  sequence_id = region_hits$candidate_id,
  sequence_type = region_hits$sequence_type,
  construct_label = region_hits$construct_label,
  round = region_hits$round,
  sequencing_batch = region_hits$sequencing_batch,
  reference_original_id = region_hits$reference_original_id,
  nucleotide_coordinates_available = region_hits$nucleotide_coordinates_available,
  protein_length_aa = region_hits$selected_protein_aa,
  gnat_region_number = region_hits$gnat_region_number,
  gnat_region_count = region_hits$gnat_region_count,
  is_primary_region = region_hits$is_primary_region,
  representative_pfam_name = region_hits$pfam_name,
  representative_pfam_accession = region_hits$pfam_accession,
  pfam_description = region_hits$pfam_description,
  aa_alignment_start = region_hits$aa_alignment_start,
  aa_alignment_end = region_hits$aa_alignment_end,
  aa_envelope_start = region_hits$aa_envelope_start,
  aa_envelope_end = region_hits$aa_envelope_end,
  domain_length_aa = region_hits$domain_length_aa,
  cds_nt_alignment_start = region_hits$cds_nt_alignment_start,
  cds_nt_alignment_end = region_hits$cds_nt_alignment_end,
  cds_nt_envelope_start = region_hits$cds_nt_envelope_start,
  cds_nt_envelope_end = region_hits$cds_nt_envelope_end,
  sequenced_insert_nt_alignment_start = region_hits$sequenced_insert_nt_alignment_start,
  sequenced_insert_nt_alignment_end = region_hits$sequenced_insert_nt_alignment_end,
  sequenced_insert_nt_envelope_start = region_hits$sequenced_insert_nt_envelope_start,
  sequenced_insert_nt_envelope_end = region_hits$sequenced_insert_nt_envelope_end,
  domain_i_evalue = region_hits$domain_i_evalue,
  domain_score = region_hits$domain_score,
  mean_posterior_probability = region_hits$mean_posterior_probability,
  overlapping_pfam_models = region_hits$overlapping_pfam_models,
  stringsAsFactors = FALSE
)

all_output <- merge(
  metadata[, c("candidate_id", "sequence_type", "construct_label", "round", "sequencing_batch", "reference_original_id", "nucleotide_coordinates_available", "used_start_1based_inclusive")],
  gnat_hits,
  by = "candidate_id",
  all.y = TRUE,
  sort = FALSE
)
all_output$cds_nt_alignment_start <- (all_output$aa_alignment_start - 1L) * 3L + 1L
all_output$cds_nt_alignment_end <- all_output$aa_alignment_end * 3L
all_output$cds_nt_alignment_start[!all_output$nucleotide_coordinates_available] <- NA_integer_
all_output$cds_nt_alignment_end[!all_output$nucleotide_coordinates_available] <- NA_integer_
all_output$sequenced_insert_nt_alignment_start <- all_output$used_start_1based_inclusive + all_output$cds_nt_alignment_start - 1L
all_output$sequenced_insert_nt_alignment_end <- all_output$used_start_1based_inclusive + all_output$cds_nt_alignment_end - 1L
all_output <- all_output[order(
  match(all_output$candidate_id, metadata$candidate_id),
  all_output$domain_i_evalue,
  -all_output$domain_score
), , drop = FALSE]

write.table(
  region_result,
  file.path(output_dir, "gnat_domain_coordinates_20.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  result,
  file.path(output_dir, "gnat_primary_domain_coordinates_20.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  all_output,
  file.path(output_dir, "gnat_all_pfam_domain_hits_20.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

md_header <- paste(
  "| Sequence | Protein length (aa) | Best Pfam GNAT model | Protein alignment (aa) | Protein envelope (aa) | Selected CDS alignment (nt) | Sequenced insert alignment (nt) | i-Evalue |",
  "|---|---:|---|---:|---:|---:|---:|---:|",
  sep = "\n"
)
format_coordinate_pair <- function(start, end) {
  if (is.na(start) || is.na(end)) "NA" else sprintf("%d–%d", start, end)
}
md_rows <- vapply(seq_len(nrow(region_result)), function(i) {
  x <- region_result[i, ]
  sprintf(
    "| %s%s | %d | %s (%s) | %s | %s | %s | %s | %.2g |",
    x$sequence_id,
    ifelse(x$gnat_region_count > 1L, paste0(" region ", x$gnat_region_number), ""),
    x$protein_length_aa, x$representative_pfam_name, x$representative_pfam_accession,
    format_coordinate_pair(x$aa_alignment_start, x$aa_alignment_end),
    format_coordinate_pair(x$aa_envelope_start, x$aa_envelope_end),
    format_coordinate_pair(x$cds_nt_alignment_start, x$cds_nt_alignment_end),
    format_coordinate_pair(x$sequenced_insert_nt_alignment_start, x$sequenced_insert_nt_alignment_end),
    x$domain_i_evalue
  )
}, character(1))

readme <- c(
  "# GNAT domain coordinates for 18 sequenced candidates plus SatA and SatB",
  "",
  "The 18 experimentally sequenced protein sequences and the SatA/SatB tree-reference proteins were scanned together against the local Pfam-A database with HMMER 3.4 `hmmscan --cut_ga`.",
  "Overlapping Pfam GNAT-family models were consolidated into physical domain regions. Nineteen proteins contain one predicted GNAT region; HIT_000035 contains two non-overlapping predicted GNAT regions.",
  "",
  "Coordinate definitions:",
  "",
  "- `aa_alignment_start/end`: residues aligned to the Pfam HMM; use these as the primary GNAT-domain coordinates.",
  "- `aa_envelope_start/end`: the wider HMMER domain envelope; useful as a conservative boundary.",
  "- `cds_nt_*`: 1-based coordinates within the selected true candidate CDS.",
  "- `sequenced_insert_nt_*`: 1-based coordinates within the original sequenced candidate insert, accounting for user-selected offsets such as nt47 or nt9.",
  "- SatA and SatB were supplied as protein-only references, so their nucleotide-coordinate fields are `NA` rather than inferred.",
  "- Coordinates are profile-based predictions, not experimentally determined structural boundaries.",
  "",
  md_header,
  md_rows,
  "",
  "Files:",
  "",
  "- `gnat_domain_coordinates_20.tsv`: all distinct GNAT regions (21 regions across 20 proteins).",
  "- `gnat_primary_domain_coordinates_20.tsv`: one strongest GNAT-domain call per protein.",
  "- `gnat_all_pfam_domain_hits_20.tsv`: every overlapping GA-threshold GNAT-family Pfam hit.",
  "- `pfam_selected_true_18_plus_SatA_SatB.domtblout`: raw HMMER domain table.",
  "- `selected_true_18_plus_SatA_SatB.faa`: the exact 20-protein analysis input.",
  "- `reference_SatA_SatB.faa`: SatA and SatB extracted from the 88-reference tree FASTA.",
  "- `sequence_metadata_20.tsv`: candidate/reference source metadata.",
  "- `run_gnat_domain_mapping.R`: reproducible analysis script."
)
writeLines(readme, file.path(output_dir, "README.md"))

summary_lines <- c(
  paste("Sequences:", nrow(result)),
  paste("Sequenced candidates:", sum(result$sequence_type == "Sequenced candidate")),
  paste("References:", sum(result$sequence_type == "Reference")),
  paste("Sequences with GA-threshold GNAT hit:", sum(!is.na(result$best_gnat_pfam_name))),
  paste("Distinct non-overlapping GNAT regions:", nrow(region_result)),
  paste("Total GA-threshold GNAT-family Pfam hits:", nrow(gnat_hits)),
  paste("Output:", file.path(output_dir, "gnat_domain_coordinates_20.tsv"))
)
writeLines(summary_lines, file.path(output_dir, "run_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
