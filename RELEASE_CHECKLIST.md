# Manuscript release checklist

- Final software name established as `LDIAG` across the package, API, CLI,
  documentation, tests, and release artifact.
- Create the public repository and add its URL to `DESCRIPTION` and `CITATION.cff`.
- Run all configured traits from a clean output directory on the cluster.
- Regenerate GWAS-pruned inputs with dosage-based reference MAF.
- Confirm section 5 reproduces the notebook's untransformed-intercept fitting
  convention on frozen study inputs.
- Compare regenerated coefficients and P-values with the manuscript tables.
- Set `plots.fail_on_error: true`; require zero `error` rows and explain every
  `skipped` row in `figures/figure_manifest.tsv`.
- Visually inspect all six manuscript PNGs and use the vector PDFs for final
  layout when the journal workflow permits.
- Render `publication/Fetal_Tissue_publication.Rmd` and archive the report.
- Run `tests/test-mixed-build-integration.R` with the release chain and add it to CI.
- Capture `sessionInfo()`, package versions, seeds, runtime, memory, and Slurm logs.
- Run `R CMD build` and `R CMD check` on the release tarball.
- Create an immutable Git tag and GitHub release.
- Archive the exact release and test data on Zenodo.
- Add the version-specific DOI to `CITATION.cff`, README, and manuscript.
