#!/usr/bin/env Rscript

if (!requireNamespace("LDIAG", quietly = TRUE)) {
  stop("Install the package first with: R CMD INSTALL .", call. = FALSE)
}
LDIAG::ldiag_cli(c("run", commandArgs(trailingOnly = TRUE)))
