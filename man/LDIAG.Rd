\name{LDIAG-package}
\alias{LDIAG-package}
\alias{LDIAG}
\alias{add_blacklist_ratio}
\alias{add_frip}
\alias{adjust_pvalues_within_trait}
\alias{build_group_accessibility}
\alias{build_peak_accessibility}
\alias{build_peak_ld_model}
\alias{build_snp_accessibility}
\alias{combine_component_tests}
\alias{compute_inverse_sqrt}
\alias{compute_snp_ld_covariance}
\alias{filter_atac_cells}
\alias{filter_atac_peaks}
\alias{fit_whitened_groups}
\alias{harmonize_gwas}
\alias{make_one_sided_wrs}
\alias{lift_genomic_ranges}
\alias{normalize_genotype}
\alias{normalize_genome_build}
\alias{read_ldiag_config}
\alias{run_atac_stage}
\alias{run_gwas_stage}
\alias{run_accessibility_stage}
\alias{run_ld_stage}
\alias{run_model_stage}
\alias{run_plot_stage}
\alias{run_ldiag}
\alias{run_wrs_stage}
\alias{ldiag_cli}
\alias{summarize_wrs_stability}
\alias{validate_ldiag_config}
\alias{validate_genome_coordinates}
\title{LDIAG: LD-aware scATAC-seq and GWAS integration}
\description{
An R workflow for ATAC quality control, GWAS harmonization, SNP- and peak-level
accessibility construction, local LD modeling, Wilcoxon scoring, LD-whitened
association testing, and aggregated Cauchy combination.
The optional plot stage writes publication figures in PNG and PDF formats.
}
\details{
See the package README and the example YAML configuration for the full input and
output contract. Block permutation and S-LDSC are outside this release.
}
\docType{package}
\keyword{package}
