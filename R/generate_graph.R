norm_vec <- function(vector, min, max) {
    vector_range <- range(vector)

    if (vector_range[1] == vector_range[2]) {
        fac <- 1
    } else {
        fac <- (max - min) / (vector_range[2] - vector_range[1])
    }

    (vector - vector_range[1]) * fac + min
}

data_paths <- function(analysis_path, actives = FALSE) {
    if (actives) {
        actives_path <- paste0(analysis_path, 'actives.csv')
    } else {
        actives_path <- NULL
    }

    matrix_data_path <- paste0(analysis_path, 'mat01.csv')
    selected_col <- paste0(analysis_path, 'selected.csv')

    list(matrix_path = matrix_data_path, selected = selected_col,
        actives_path = actives_path)
}

load_data <- function(parameters, analysis_path) {
    if (parameters$type != 'simimatrix' ||
        parameters$type == 'simiclustermatrix') {

        paths_list <- data_paths(analysis_path, actives = TRUE)

        matrix_data <- Matrix::readMM(paths_list$matrix_path)
        actives <- read.table(paths_list$actives_path, sep = '\t', quote = '"')

        colnames(matrix_data) <- actives[, 1]
    } else if (parameters$type == 'simimatrix' ||
        parameters$type == 'simiclustermatrix') {

        paths_list <- data_paths(analysis_path)

        matrix_data <- read.csv2(paths_list$matrix_path)
        matrix_data <- as.matrix(matrix_data)
    }

    list(paths = paths_list, matrix = matrix_data)
}

load_coords <- function() {
    coords <- try(coords, TRUE)
    load(paste0(analysis_path, 'RData.RData'))

    if (!is.matrix(coords)) {
        coords <- NULL
    }

    coords
}

select_matrix_word <- function(parameters, matrix_data, selected_column) {
    if ('word' %in% parameters) {
        word <- TRUE
        index <- parameters$word + 1
    } else {
        word <- FALSE
        index <- NULL
    }

    if (!word) {
        matrix_data <- matrix_data[, selected_column]
    } else {
        forme <- colnames(matrix_data)[index]

        if (!index %in% selected_column) {
            selected_column <- append(selected_column, index)
        }

        matrix_data <- matrix_data[, selected_column]
        index <- which(colnames(matrix_data) == forme)
    }

    list(matrix = matrix_data, index = index)
}

check_selected_column <- function(sel_col_file, matrix_data) {
    if (file.exists(sel_col_file)) {
        sel_col <- read.csv2(sel_col_file, header = FALSE)
        sel_col <- sel_col[, 1] + 1
    } else {
        sel_col <- 1:ncol(matrix_data)
    }

    sel_col
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

check_word_parameter <- function(parameters, sparse, matrix_data, index) {
    if ('word' %in% parameters) {
        sparse <- graph_word(sparse, index)
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
generate_graph <- function(parameters, analysis_path) {
    parameters <- check_coeff_tv(parameters)
    parameters <- check_vcex(parameters)

    if (parameters$keep_coord) {
        parameters$coords <- load_coords()
    } else {
        parameters$coords <- NULL
    }

    data_definition <- load_data(parameters, analysis_path)

    selected_colmun <- check_selected_column(data_definition$paths$selected,
        data_definition$matrix)
    selection <- select_matrix_word(parameters, data_definition$matrix, selected_colmun)

    sparse <- create_sparse(parameters, selection$matrix)
    sparse <- check_inf(sparse)

    valid_data <- check_word_parameter(parameters, sparse, selection$matrix,
        selection$index)

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
