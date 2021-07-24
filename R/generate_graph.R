norm_vec <- function(vector, min, max) {
    vector_range <- range(vector)

    if (vector_range[1] == vector_range[2]) {
        fac <- 1
    } else {
        fac <- (max - min) / (vector_range[2] - vector_range[1])
    }

    (vector - vector_range[1]) * fac + min
}

load_coords <- function() {
    coords <- try(coords, TRUE)
    load(paste0(analysis_path, 'RData.RData'))

    if (!is.matrix(coords)) {
        coords <- NULL
    }

    coords
}

check_inf <- function(sparse) {
    if (length(which(sparse == Inf))) {
        infp <- which(sparse == Inf)
        sparse[infp] <- NA
        maxmat <- max(sparse, na.rm = TRUE)

        if (maxmat > 0) {
            maxmat <- maxmat + 1
        } else {
            maxmat <- 0
        }

        sparse[infp] <- maxmat
    }

    if (length(which(sparse == -Inf))) {
        infm <- which(sparse == -Inf)
        sparse[infm] <- NA
        minmat <- min(sparse, na.rm = TRUE)

        if (maxmat < 0) {
            minmat <- minmat - 1
        } else {
            minmat <- 0
        }

        sparse[infm] <- minmat
    }

    sparse
}

check_word_parameter <- function(parameters, sparse, matrix_data) {
    if ('word' %in% parameters) {
        sparse <- graph_word(sparse, parameters$word)
        col_sum <- colSums(sparse)

        if (length(col_sum)) sparse <- sparse[, -which(col_sum == 0)]

        row_sum <- rowSums(sparse)

        if (length(row_sum)) sparse <- sparse[-which(row_sum == 0),]
        if (length(col_sum)) matrix_data <- matrix_data[, -which(col_sum == 0)]
    }

    list(sparse = sparse, matrix = matrix_data)
}

create_sparse <- function(parameters, matrix_data) {
    if (parameters$method == 'cooc') {
        sparse <- square_matrix(matrix_data)

    } else if (parameters$method == 'Russel') {
        sparse <- proxy::simil(matrix_data, method = parameters$method,
            diag = TRUE, upper = TRUE, by_rows = FALSE)

    } else if (parameters$method == 'binomial') {
        sparse <- binom_sim(matrix_data)

    } else {
        sparse <- proxy::simil(as.matrix(matrix_data), method = parameters$method,
            diag = TRUE, upper = TRUE, by_rows = FALSE)
    }

    as.matrix(stats::as.dist(sparse, diag = TRUE, upper = TRUE))
}

check_coeff_tv <- function(parameters) {
    if (parameters$coeff_tv && parameters$sfromchi) {
        parameters$coeff_tv <- NULL
        parameters$minmaxeff <- c(parameters$tvmin, parameters$tvmax)
    } else {
        parameters$minmaxeff <- NULL
    }

    parameters
}

check_vcex <- function(parameters) {
    if (parameters$vcex || parameters$cexfromchi) {
        parameters$vcexminmax <- c(parameters$vcexmin, parameters$vcexmax)
    } else {
        parameters$vcexminmax <- NULL
    }

    parameters
}

#' Generate similitude graph
#' 
#' @param parameters list
#' @param analysis_path of data
#' @return list of objects
#' 
#' @export
generate_graph <- function(parameters, dtm) {
    parameters <- check_coeff_tv(parameters)
    parameters <- check_vcex(parameters)

    if (parameters$keep_coord) {
        parameters$coords <- load_coords()
    } else {
        parameters$coords <- NULL
    }

    sparse <- create_sparse(parameters, dtm)
    sparse <- check_inf(sparse)

    valid_data <- check_word_parameter(parameters, sparse, dtm)

    eff <- colSums(as.matrix(valid_data$matrix))
    x <- list(mat = valid_data$sparse, eff = eff)

    do_simi(
        x,
        method = parameters$method,
        seuil = parameters$seuil,
        plot_type = parameters$type,
        layout_type = parameters$layout,
        max_tree = parameters$max_tree,
        coeff_vertex = parameters$coeff_tv,
        coeff_edge = parameters$coeff_te_range,
        minmax_eff = parameters$minmaxeff,
        vcex_minmax = parameters$vcexminmax,
        cex = parameters$cex,
        coords = parameters$coords,
        communities = parameters$communities,
        halo = parameters$halo
    )
}
