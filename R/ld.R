scale_genotype_for_ld <- function(genotype) {
  validate_named_matrix(genotype, "genotype")
  x <- t(as.matrix(genotype))
  if (nrow(x) < 2L) stop("At least two reference samples are required for LD.", call. = FALSE)
  means <- colMeans(x, na.rm = TRUE)
  missing <- which(is.na(x), arr.ind = TRUE)
  if (nrow(missing)) {
    x[missing] <- means[missing[, "col"]]
  }
  centered <- sweep(x, 2L, means, "-")
  standard_deviation <- sqrt(colSums(centered^2) / (nrow(centered) - 1))
  keep <- is.finite(standard_deviation) & standard_deviation > 0
  if (!any(keep)) stop("Every reference variant has zero or invalid variance.", call. = FALSE)
  x <- sweep(centered[, keep, drop = FALSE], 2L, standard_deviation[keep], "/")
  colnames(x) <- rownames(genotype)[keep]
  x
}

as_general_csparse <- function(matrix) {
  methods::as(methods::as(matrix, "generalMatrix"), "CsparseMatrix")
}

ld_from_scaled <- function(scaled_genotype, snps1, snps2) {
  result <- crossprod(
    scaled_genotype[, snps1, drop = FALSE],
    scaled_genotype[, snps2, drop = FALSE]
  ) / (nrow(scaled_genotype) - 1)
  result <- as.matrix(result)
  result[!is.finite(result)] <- 0
  rownames(result) <- snps1
  colnames(result) <- snps2
  result
}

#' Construct local SNP-level LD response covariance
#'
#' @param genotype Variant-by-sample dosage matrix named by genomic coordinate.
#' @param gwas Harmonized GWAS table containing `coord_id`, chromosome, and position.
#' @param columns GWAS column mapping.
#' @param window_bp Maximum same-chromosome distance represented in the covariance.
#' @param factor_two Use `2 * R^2` instead of the response correlation `R^2`.
#' @param correlation_block_size Number of focal variants handled in each BLAS cross-product.
#' @return A sparse symmetric covariance matrix ordered by chromosome and position.
#' @export
compute_snp_ld_covariance <- function(genotype, gwas,
                                      columns = ldiag_defaults()$gwas$columns,
                                      window_bp = 5000000L, factor_two = FALSE,
                                      correlation_block_size = 256L) {
  require_namespaces("Matrix", "SNP LD covariance")
  if (!"coord_id" %in% colnames(gwas)) {
    gwas$coord_id <- coordinate_id(gwas[[columns$chr]], gwas[[columns$pos]])
  }
  order_ids <- ordered_intersection(gwas$coord_id, rownames(genotype))
  gwas <- gwas[match(order_ids, gwas$coord_id), , drop = FALSE]
  genotype <- genotype[order_ids, , drop = FALSE]
  ord <- order(chr_to_order(gwas[[columns$chr]]), as.integer(gwas[[columns$pos]]))
  gwas <- gwas[ord, , drop = FALSE]
  genotype <- genotype[ord, , drop = FALSE]
  scaled <- scale_genotype_for_ld(genotype)
  keep <- gwas$coord_id %in% colnames(scaled)
  gwas <- gwas[keep, , drop = FALSE]
  scaled <- scaled[, gwas$coord_id, drop = FALSE]
  n <- nrow(gwas)
  if (n == 0L) stop("No variable reference variants remain for SNP LD.", call. = FALSE)

  i_list <- vector("list", n)
  j_list <- vector("list", n)
  x_list <- vector("list", n)
  counter <- 0L
  by_chr <- split(seq_len(n), chr_to_order(gwas[[columns$chr]]))
  for (index in by_chr) {
    position <- as.integer(gwas[[columns$pos]][index])
    ids <- gwas$coord_id[index]
    x_chr <- scaled[, ids, drop = FALSE]
    block_starts <- seq.int(1L, length(index), by = as.integer(correlation_block_size))
    for (block_start in block_starts) {
      block_end <- min(length(index), block_start + as.integer(correlation_block_size) - 1L)
      focal <- seq.int(block_start, block_end)
      maximum_right <- findInterval(position[[block_end]] + window_bp, position)
      comparison <- seq.int(block_start, maximum_right)
      correlation_block <- crossprod(
        x_chr[, focal, drop = FALSE], x_chr[, comparison, drop = FALSE]
      ) / (nrow(x_chr) - 1)
      for (offset in seq_along(focal)) {
        a <- focal[[offset]]
        right <- findInterval(position[[a]] + window_bp, position)
        js <- seq.int(a, right)
        block_columns <- seq.int(a - block_start + 1L, right - block_start + 1L)
        value <- as.numeric(correlation_block[offset, block_columns, drop = TRUE])^2
        if (factor_two) value <- 2 * value
        counter <- counter + 1L
        i_list[[counter]] <- rep.int(index[[a]], length(js))
        j_list[[counter]] <- index[js]
        x_list[[counter]] <- value
      }
    }
  }
  sigma <- Matrix::sparseMatrix(
    i = unlist(i_list[seq_len(counter)], use.names = FALSE),
    j = unlist(j_list[seq_len(counter)], use.names = FALSE),
    x = unlist(x_list[seq_len(counter)], use.names = FALSE),
    dims = c(n, n),
    dimnames = list(gwas$coord_id, gwas$coord_id)
  )
  sigma <- sigma + Matrix::t(sigma) - Matrix::Diagonal(x = Matrix::diag(sigma))
  Matrix::diag(sigma) <- if (factor_two) 2 else 1
  as_general_csparse(sigma)
}

#' Compute a ridge-stabilized inverse square root
#'
#' @param matrix Symmetric square matrix.
#' @param ridge Nonnegative ridge added to every eigenvalue.
#' @param eigen_tolerance Lower eigenvalue floor.
#' @return Dense symmetric inverse square-root matrix.
#' @export
compute_inverse_sqrt <- function(matrix, ridge = 0.01, eigen_tolerance = 1e-8) {
  matrix <- as.matrix(matrix)
  if (nrow(matrix) == 0L || nrow(matrix) != ncol(matrix)) {
    stop("matrix must be a nonempty square matrix.", call. = FALSE)
  }
  matrix <- (matrix + t(matrix)) / 2
  eig <- eigen(matrix, symmetric = TRUE)
  value <- pmax(eig$values, eigen_tolerance) + ridge
  result <- eig$vectors %*% diag(1 / sqrt(value), nrow = length(value)) %*% t(eig$vectors)
  rownames(result) <- rownames(matrix)
  colnames(result) <- colnames(matrix)
  result
}

split_sparse_components <- function(matrix) {
  require_namespaces(c("Matrix", "igraph"), "LD connected components")
  entries <- Matrix::summary(matrix)
  entries <- entries[entries$i != entries$j & entries$x != 0, , drop = FALSE]
  nodes <- rownames(matrix)
  if (nrow(entries) == 0L) return(setNames(as.list(nodes), paste0("block_", seq_along(nodes))))
  edges <- cbind(nodes[entries$i], nodes[entries$j])
  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  missing <- setdiff(nodes, igraph::V(graph)$name)
  if (length(missing)) graph <- igraph::add_vertices(graph, nv = length(missing), name = missing)
  membership <- igraph::components(graph)$membership
  groups <- split(names(membership), membership)
  groups <- lapply(groups, function(group) nodes[nodes %in% group])
  names(groups) <- paste0("block_", seq_along(groups))
  groups
}

inverse_sqrt_components <- function(sigma, ridge, eigen_tolerance) {
  require_namespaces("Matrix", "block inverse square roots")
  components <- split_sparse_components(sigma)
  matrices <- lapply(components, function(ids) {
    compute_inverse_sqrt(sigma[ids, ids, drop = FALSE], ridge, eigen_tolerance)
  })
  combined <- Matrix::bdiag(matrices)
  combined <- as_general_csparse(combined)
  ids <- unlist(lapply(matrices, rownames), use.names = FALSE)
  rownames(combined) <- ids
  colnames(combined) <- ids
  list(matrix = combined, blocks = matrices, indices = components)
}

grow_numeric_storage <- function(i, j, x, needed) {
  if (needed <= length(i)) return(list(i = i, j = j, x = x))
  new_length <- max(needed, max(1024L, length(i) * 2L))
  length(i) <- new_length
  length(j) <- new_length
  length(x) <- new_length
  list(i = i, j = j, x = x)
}

#' Build peak-level GWAS response and local LD covariance
#'
#' @param gwas Harmonized GWAS table containing `coord_id` and `Z`.
#' @param genotype Variant-by-sample genotype matrix named by `coord_id`.
#' @param peak_matrix Peak-by-cell accessibility matrix.
#' @param overlap Data frame with `snp` and `peak` columns.
#' @param window_bp Local peak-pair distance.
#' @param ridge Ridge used for within-peak and response covariance matrices.
#' @param eigen_tolerance Eigenvalue floor.
#' @param max_snps_for_chr_ld Maximum variants for caching a chromosome LD matrix.
#' @return List containing peak response, covariance, and inverse square root.
#' @export
build_peak_ld_model <- function(gwas, genotype, peak_matrix, overlap,
                                window_bp = 5000000L, ridge = 0.01,
                                eigen_tolerance = 1e-8,
                                max_snps_for_chr_ld = 15000L) {
  require_namespaces(c("Matrix", "Signac", "GenomicRanges"), "peak LD model")
  validate_named_matrix(genotype, "peak-model genotype")
  validate_named_matrix(peak_matrix, "peak accessibility")
  if (!all(c("snp", "peak") %in% colnames(overlap))) {
    stop("overlap must contain snp and peak columns.", call. = FALSE)
  }
  if (!all(c("coord_id", "Z") %in% colnames(gwas))) {
    stop("gwas must contain coord_id and Z columns.", call. = FALSE)
  }
  scaled <- scale_genotype_for_ld(genotype)
  z_lookup <- setNames(as.numeric(gwas$Z), gwas$coord_id)
  overlap <- unique(overlap[
    overlap$snp %in% colnames(scaled) & overlap$snp %in% names(z_lookup) &
      overlap$peak %in% rownames(peak_matrix),
    c("snp", "peak"), drop = FALSE
  ])
  if (nrow(overlap) == 0L) stop("No valid SNP-peak overlaps remain for peak LD.", call. = FALSE)
  peak_snps <- lapply(split(overlap$snp, overlap$peak), unique)

  peak_info <- lapply(peak_snps, function(snps) {
    snps <- intersect(snps, colnames(scaled))
    z <- z_lookup[snps]
    keep <- is.finite(z)
    snps <- snps[keep]
    z <- as.numeric(z[keep])
    if (length(snps) == 0L) return(NULL)
    r <- ld_from_scaled(scaled, snps, snps)
    inverse <- compute_inverse_sqrt(r, ridge, eigen_tolerance)
    y <- sum(as.numeric(inverse %*% z)^2) / length(snps)
    if (!is.finite(y)) return(NULL)
    list(snps = snps, k = length(snps), y = y, inverse = inverse)
  })
  peak_info <- peak_info[!vapply(peak_info, is.null, logical(1))]
  peak_ids <- rownames(peak_matrix)[rownames(peak_matrix) %in% names(peak_info)]
  peak_info <- peak_info[peak_ids]
  y <- vapply(peak_info, `[[`, numeric(1), "y")
  if (length(y) == 0L) stop("No valid peak response values could be computed.", call. = FALSE)

  peak_ranges <- Signac::StringToGRanges(peak_ids, sep = c("-", "-"))
  names(peak_ranges) <- peak_ids
  peak_table <- data.frame(
    peak = peak_ids,
    chr = as.character(GenomicRanges::seqnames(peak_ranges)),
    start = GenomicRanges::start(peak_ranges),
    end = GenomicRanges::end(peak_ranges),
    stringsAsFactors = FALSE
  )
  peak_table$chr_order <- chr_to_order(peak_table$chr)
  peak_table <- peak_table[order(peak_table$chr_order, peak_table$start, peak_table$end), ]
  peak_index <- setNames(seq_along(peak_ids), peak_ids)

  capacity <- max(1024L, length(peak_ids) * 4L)
  pair_i <- integer(capacity)
  pair_j <- integer(capacity)
  pair_x <- numeric(capacity)
  count <- 0L
  for (chromosome in sort(unique(peak_table$chr_order))) {
    table_chr <- peak_table[peak_table$chr_order == chromosome, , drop = FALSE]
    chromosome_snps <- unique(unlist(lapply(table_chr$peak, function(p) peak_info[[p]]$snps), FALSE))
    chromosome_snps <- intersect(chromosome_snps, colnames(scaled))
    cached <- if (length(chromosome_snps) <= max_snps_for_chr_ld) {
      ld_from_scaled(scaled, chromosome_snps, chromosome_snps)
    } else NULL
    right_boundary <- findInterval(table_chr$end + window_bp, table_chr$start)

    for (a in seq_len(nrow(table_chr))) {
      if (right_boundary[[a]] < a) next
      p1 <- table_chr$peak[[a]]
      info1 <- peak_info[[p1]]
      for (b in seq.int(a, right_boundary[[a]])) {
        p2 <- table_chr$peak[[b]]
        info2 <- peak_info[[p2]]
        cross_ld <- if (is.null(cached)) {
          ld_from_scaled(scaled, info1$snps, info2$snps)
        } else cached[info1$snps, info2$snps, drop = FALSE]
        transformed <- info1$inverse %*% cross_ld %*% info2$inverse
        covariance_q <- 2 * sum(transformed * transformed)
        if (!is.finite(covariance_q) || covariance_q == 0) next
        count <- count + 1L
        storage <- grow_numeric_storage(pair_i, pair_j, pair_x, count)
        pair_i <- storage$i
        pair_j <- storage$j
        pair_x <- storage$x
        pair_i[[count]] <- peak_index[[p1]]
        pair_j[[count]] <- peak_index[[p2]]
        pair_x[[count]] <- covariance_q
      }
    }
  }
  if (count == 0L) stop("No peak-pair covariance entries were generated.", call. = FALSE)
  covariance_q <- Matrix::sparseMatrix(
    i = pair_i[seq_len(count)], j = pair_j[seq_len(count)], x = pair_x[seq_len(count)],
    dims = c(length(peak_ids), length(peak_ids)), dimnames = list(peak_ids, peak_ids)
  )
  covariance_q <- covariance_q + Matrix::t(covariance_q) - Matrix::Diagonal(x = Matrix::diag(covariance_q))
  k <- vapply(peak_info, `[[`, numeric(1), "k")[peak_ids]
  sigma <- methods::as(methods::as(covariance_q, "generalMatrix"), "TsparseMatrix")
  denominator <- k[sigma@i + 1L] * k[sigma@j + 1L]
  sigma@x <- ifelse(is.finite(denominator) & denominator > 0, sigma@x / denominator, 0)
  sigma <- as_general_csparse(sigma)
  inverse <- inverse_sqrt_components(sigma, ridge, eigen_tolerance)
  list(
    response = y[peak_ids],
    sigma = sigma,
    covariance_q = covariance_q,
    inverse_sqrt = inverse$matrix,
    inverse_sqrt_blocks = inverse$blocks,
    block_indices = inverse$indices,
    peak_info = lapply(peak_info, function(info) info[c("snps", "k")])
  )
}

align_genotype_to_coordinates <- function(gwas, genotype, columns) {
  index <- match(as.character(gwas[[columns$rsid]]), rownames(genotype))
  if (anyNA(index)) stop("Harmonized genotype does not align to GWAS rsIDs.", call. = FALSE)
  genotype <- genotype[index, , drop = FALSE]
  rownames(genotype) <- gwas$coord_id
  genotype
}

#' Build local SNP- and peak-level LD models
#'
#' @param config Parsed LDIAG configuration.
#' @param traits Optional trait subset.
#' @return Named list of generated files.
#' @export
run_ld_stage <- function(config, traits = names(config$gwas$traits)) {
  output_dir <- stage_output_dir(config, "03_ld")
  outputs <- setNames(vector("list", length(traits)), traits)
  for (trait in traits) {
    columns <- deep_merge(config$gwas$columns, trait_config(config, trait)$columns %||% list())
    gwas <- load_serialized(file.path(config$output_dir, "02_gwas", paste0(trait, ".gwas.rds")))
    genotype <- load_serialized(file.path(config$output_dir, "02_gwas", paste0(trait, ".genotype.rds")))
    genotype <- align_genotype_to_coordinates(gwas, genotype, columns)
    snp_matrix <- load_serialized(file.path(config$output_dir, "02_accessibility", paste0(trait, ".snp_by_cell.rds")))
    peak_matrix <- load_serialized(file.path(config$output_dir, "02_accessibility", paste0(trait, ".peak_by_cell.rds")))
    overlap <- utils::read.delim(
      file.path(config$output_dir, "02_accessibility", paste0(trait, ".snp_peak_overlap.tsv")),
      stringsAsFactors = FALSE, check.names = FALSE
    )

    snp_ids <- ordered_intersection(rownames(snp_matrix), gwas$coord_id, rownames(genotype))
    gwas_snp <- gwas[match(snp_ids, gwas$coord_id), , drop = FALSE]
    genotype_snp <- genotype[snp_ids, , drop = FALSE]
    snp_sigma <- compute_snp_ld_covariance(
      genotype_snp, gwas_snp, columns, config$ld$window_bp,
      correlation_block_size = config$ld$correlation_block_size
    )
    snp_inverse <- inverse_sqrt_components(
      snp_sigma, config$ld$ridge, config$ld$eigen_tolerance
    )
    peak_model <- build_peak_ld_model(
      gwas, genotype, peak_matrix, overlap,
      window_bp = config$ld$window_bp,
      ridge = config$ld$ridge,
      eigen_tolerance = config$ld$eigen_tolerance,
      max_snps_for_chr_ld = config$ld$max_snps_for_chr_ld %||% 15000L
    )

    paths <- list(
      snp_sigma = file.path(output_dir, paste0(trait, ".snp_sigma.rds")),
      snp_inverse = file.path(output_dir, paste0(trait, ".snp_inverse_sqrt.rds")),
      snp_blocks = file.path(output_dir, paste0(trait, ".snp_blocks.rds")),
      peak_model = file.path(output_dir, paste0(trait, ".peak_ld_model.rds"))
    )
    save_rds_atomic(snp_sigma, paths$snp_sigma)
    save_rds_atomic(snp_inverse$matrix, paths$snp_inverse)
    save_rds_atomic(snp_inverse$indices, paths$snp_blocks)
    save_rds_atomic(peak_model, paths$peak_model)
    outputs[[trait]] <- paths
    message("LD stage complete for ", trait)
  }
  outputs
}
