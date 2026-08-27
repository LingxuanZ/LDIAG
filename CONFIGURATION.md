# Configuration reference

All paths may be absolute or relative to the YAML configuration file. Values
shown below are defaults unless marked as required.

## Project and execution

| Key | Default | Purpose |
| --- | --- | --- |
| `project.name` | `LDIAG analysis` | Human-readable run name. |
| `project.output_dir` | `results` | Root for all generated files. |
| `parameters.seed` | `1000` | Reproducible LD pruning and balanced sampling seed. |
| `parameters.workers` | `1` | BiocParallel workers for fragment overlap calculations. |

## Genome harmonization

| Key | Default | Purpose |
| --- | --- | --- |
| `genome.target` | `hg19` | Common build for saved GWAS, SNP, peak, LD, and model inputs; accepts hg19/GRCh37 or hg38/GRCh38. |
| `genome.chains.hg19_to_hg38` | null | Directional UCSC hg19-to-hg38 chain file. |
| `genome.chains.hg38_to_hg19` | null | Directional UCSC hg38-to-hg19 chain file. |
| `genome.unmapped_action` | `drop` | `drop` or `error` for ranges with no target mapping. |
| `genome.multimapping_action` | `drop` | `drop` or `error` for ranges with multiple target mappings. |
| `genome.duplicate_action` | `drop` | `drop` or `error` when multiple inputs produce the same target range. |

Input builds are explicit because coordinates alone cannot reliably distinguish
hg19 from hg38. The software checks chromosome bounds and compares declared
builds with available Seqinfo or GWAS metadata. A conversion requires the chain
in the correct direction. ATAC-to-target is required when those builds differ.
Target-to-ATAC is additionally required only for a GWAS that did not originate
in the ATAC build; otherwise its preserved input coordinates are reused for
source-build fragment queries.

The standard UCSC files are `hg19ToHg38.over.chain.gz` and
`hg38ToHg19.over.chain.gz`. Keep local immutable copies and record their MD5;
LDIAG also writes the observed chain MD5 into its QC output.

## ATAC input and QC

| Key | Default | Purpose |
| --- | --- | --- |
| `atac.input` | required | Seurat/Signac `.rds` or `.rda` object. |
| `atac.object_name` | null | Object name when an `.rda` contains multiple objects. |
| `atac.assay` | `ATAC` | Chromatin assay name. |
| `atac.genome` | `hg19` | Declared build of both ATAC peaks and fragment files; hg19/GRCh37 or hg38/GRCh38. |
| `atac.tissue_column` | `tissue` | Tissue metadata column. |
| `atac.group_column` | `tissue_cell_type` | Column used for minimum group-size QC. |
| `atac.attach_fragments` | `false` | Attach files listed under `atac.fragment_files`. |
| `atac.fragment_files` | null | Named tissue-to-fragment-file mapping. |
| `atac.validate_fragments` | `true` | Ask Signac to validate cell barcodes during attachment. |
| `atac.compute_qc_metrics` | `false` | Compute nucleosome, TSS, FRIP, and blacklist metrics. |
| `atac.recompute_embedding` | `false` | Re-run TF-IDF, LSI, UMAP, neighbors, and clustering after QC. |
| `atac.qc.tss_min` | `2.5` | Minimum TSS enrichment. |
| `atac.qc.frip_min` | `0.3` | Minimum fraction of fragments in peaks. |
| `atac.qc.count_quantile_min` | `0.01` | Remove cells below this ATAC-count quantile. |
| `atac.qc.blacklist_max` | `0.02` | Exclusive upper blacklist-fragment ratio. |
| `atac.qc.nucleosome_max` | `4` | Maximum nucleosome signal. |
| `atac.qc.min_cells_per_group` | `50` | Minimum retained cells per annotation group. |
| `atac.qc.peak_width_min` | `100` | Minimum peak width in base pairs. |
| `atac.qc.peak_width_max` | `1000` | Maximum peak width in base pairs. |
| `atac.qc.peak_min_total_counts` | `1` | Minimum total counts across cells for a peak. |
| `atac.qc.peak_min_cells` | `1` | Minimum number of cells in which a peak is detected. |
| `atac.qc.peak_autosomes_only` | `true` | Retain peaks on chromosomes 1-22 only. |
| `atac.qc.peak_remove_blacklist` | `true` | Remove peaks using the Signac blacklist matching `atac.genome`. |

When `compute_qc_metrics` is false, the five QC metadata fields must already be
present. Fragment objects must be attached before the accessibility stage even
if QC metrics are precomputed. The ATAC stage writes both the configured values
and the data-derived ATAC count cutoff to `01_atac/qc_thresholds.tsv`.

## GWAS and reference genotypes

| Key | Default | Purpose |
| --- | --- | --- |
| `gwas.traits.<name>.summary` | required | Aligned GWAS table for one trait. |
| `gwas.traits.<name>.genotype` | required | Matrix-like genotype object or probability array/list. |
| `summary_object` / `genotype_object` | null | Object names for multi-object `.rda` inputs. |
| `gwas.genome` | `hg19` | Default declared GWAS input build; each named GWAS may override it. |
| `gwas.genome_column` | null | Optional summary-statistics column containing build metadata to verify. |
| `gwas.maf_min` | `0.05` | Minimum MAF in both GWAS and reference panel. |
| `gwas.genotype_format` | `auto` | `dosage`, `hard_call`, `gt`, `allele`, `gp`, `gl`, `pl`, or automatic detection. |
| `gwas.genotype_orientation` | `auto` | `variants_by_samples`, `samples_by_variants`, or infer from rsIDs. |
| `gwas.genotype_ploidy` | `2` | Allele copy number used to derive reference frequency. |
| `gwas.genotype_missing_values` | `[-9]` | Additional input values decoded as missing. |
| `gwas.genotype_variant_id_column` | null | Variant ID column for genotype data-frame input. |
| `gwas.reference_missing_rate_max` | `0.01` | Maximum missing reference genotypes per variant. |
| `gwas.remove_ambiguous` | `true` | Remove A/T, T/A, C/G, and G/C variants. |
| `gwas.ld_prune` | `true` | Perform SNPRelate pruning. |
| `gwas.ld_prune_tool` | `snprelate` | LD-pruning backend; only SNPRelate is integrated in this release. |
| `gwas.ld_prune_noninteger` | `error` | Reject continuous dosage during pruning, or `round` to opt into nearest hard calls. |
| `gwas.ld_method` | `composite` | SNPRelate statistic: `composite`, `r`, `dprime`, or `corr`. |
| `gwas.ld_threshold` | `0.8` | SNPRelate LD threshold. |
| `gwas.ld_window_bp` | `500000` | LD pruning window in base pairs. |
| `gwas.ld_slide_max_n` | null | Optional maximum variants in each pruning window. |
| `gwas.ld_prune_maf_min` | `0.005` | Additional SNPRelate MAF filter during pruning. |
| `gwas.ld_remove_monomorphic` | `true` | Ask SNPRelate to remove monomorphic variants. |
| `gwas.ld_start_pos` | `random.f500` | Window start: `random.f500`, `random`, `first`, or `last`. |
| `gwas.ld_num_threads` | `1` | Threads used by `snpgdsLDpruning()`. |

`gwas.columns` maps `rsid`, `chr`, `pos`, `a1`, `a2`, `maf`, `p`, `beta`,
and `se` to input columns. A GWAS may override the column mapping and every
GWAS QC, genotype, and LD-pruning setting under `gwas.traits.<name>`. The
effective settings and variant counts are written to each
`02_gwas/<name>.gwas_qc.tsv` and to the combined QC summary.

GWAS and genotype inputs may be `.rds`, `.rda/.RData`, `.csv[.gz]`,
`.tsv[.gz]`, or `.txt[.gz]`. A text genotype table should place stable variant
IDs in the configured `genotype_variant_id_column` and samples in the remaining
columns. Genotype matrices are joined by those IDs; liftOver changes GWAS
coordinates, not dosage values. Native PLINK/VCF/BGEN/GDS containers still need
to be read or exported to a supported matrix first.

`auto` accepts numeric allele dosage, biallelic GT strings, allele labels, and
well-formed GP probabilities. GL and PL must be declared explicitly because
their comma-delimited representation cannot be distinguished safely from GP.
Logical or declared `carrier` matrices are rejected: carrier status does not
identify heterozygous versus homozygous carriers and therefore cannot recover
exact allele frequency. Numeric `0/1` with diploid ploidy is likewise ambiguous
under `auto`; declare `dosage` only for true allele counts, or set ploidy to `1`
for haploid calls.

SNPRelate pruning operates on diploid `0/1/2` hard calls. Haploid dosage is
scaled exactly to `0/2`. Continuous dosage, GP/GL/PL-derived expected dosage, and
polyploid values that do not map exactly require either
`ld_prune_noninteger: round` or `ld_prune: false`. The normalized dosage is
retained unchanged for downstream LD calculations even when pruning uses the
explicitly requested nearest hard calls.

The first argument of `snpgdsLDpruning()` is an open GDS file handle created
internally; it is not a token. The user-facing selector is `ld_prune_tool`, and
the LD statistic is selected with `ld_method`. To prune with PLINK or another
external program in this release, prepare matching pre-pruned GWAS and genotype
inputs first and set `ld_prune: false`; arbitrary external commands are not run
from YAML.

The package does not infer allele flips from a dosage matrix without allele
metadata. Summary statistics must be aligned to the reference panel before
input.

## Accessibility, LD, and tests

| Key | Default | Purpose |
| --- | --- | --- |
| `accessibility.chunks` | `20` | Chunks used for SNP FeatureMatrix calculation. |
| `ld.window_bp` | `5000000` | Same-chromosome local covariance distance. |
| `ld.ridge` | `0.01` | Ridge added before inverse square roots. |
| `ld.eigen_tolerance` | `1e-8` | Eigenvalue floor before ridge regularization. |
| `ld.correlation_block_size` | `256` | Focal variants per vectorized BLAS correlation block; performance-only setting. |
| `ld.max_snps_for_chr_ld` | `15000` | Cache full chromosome LD below this size; otherwise calculate peak-pair LD on demand. |
| `wrs.truncation_quantiles` | `[0.02, 0.98]` | Winsorization limits for enrichment-oriented Wilcoxon Z-scores. |
| `wrs.balanced.enabled` | `false` | Run repeated equal-size one-versus-rest sampling. |
| `wrs.balanced.cells_per_side` | `[500, 1000]` | Target and background sample sizes. |
| `wrs.balanced.repeats` | `30` | Repeats at each sample size. |
| `wrs.balanced.auc_cutoff` | `0.55` | Selection threshold for repeat-pair stability. |
| `wrs.balanced.stability_frequency_cutoff` | `0.8` | Minimum selected-repeat proportion for a feature to be marked stable. |
| `wrs.balanced.top_n` | `100` | Top-ranked features compared between repeat pairs. |
| `model.groupings` | `[tissue_cell_type, tissue]` | Metadata resolutions modeled independently. |
| `model.require_all_components` | `true` | Require all four component P-values for ACAT. |

Changing `ld.ridge`, `ld.window_bp`, Wilcoxon truncation, or group definitions
changes the statistical model and should be treated as a new analysis rather
than an execution-only option.

`ld.correlation_block_size` changes batching only. Block size 1 and larger
blocks use the same standardized dosage values, local-distance mask, and output
schema; numerical-equivalence tests cover this contract.

## Publication figures

| Key | Default | Purpose |
| --- | --- | --- |
| `plots.formats` | `[png, pdf]` | One or more of `png` and `pdf`. |
| `plots.dpi` | `300` | Raster resolution for PNG output. |
| `plots.base_size` | `11` | Base ggplot text size. |
| `plots.max_umap_cells` | `100000` | Deterministic maximum nuclei per UMAP panel. |
| `plots.count_bin_max` | `4` | Last exact accessibility-count bin; larger values form one `N+` bin. |
| `plots.ld_heatmap_max_features` | `250` | Maximum SNPs/peaks shown in a dense diagnostic LD heatmap. |
| `plots.wrs_histogram_bins` | `50` | WRS histogram bins. |
| `plots.top_features_per_group` | `3` | WRS-ranked features selected per trait and annotation group. |
| `plots.pvalue_cap` | `12` | Display cap for `-log10(P)`; raw table values are unchanged. |
| `plots.sldsc_summary` | null | Optional merged S-LDSC TSV/RDS/RData used for method-comparison figures. |
| `plots.fail_on_error` | `false` | Stop immediately on a plotting error; use `true` for release runs. |

The `plot` stage reads saved results and does not modify analytical objects.
Every attempted figure is recorded in `figures/figure_manifest.tsv` with its
source path, output paths, status, message, and elapsed time. S-LDSC computation
is external; accepted summaries need `trait`, an annotation/group column, and a
coefficient z-score. Common coefficient/enrichment P-value column names are
standardized automatically.
