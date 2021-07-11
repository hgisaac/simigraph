norm.vec <- function(vector, min, max) {
  vector.range <- range(vector)

  if (vector.range[1] == vector.range[2]) {
    fac <- 1
  } else {
    fac <- (max - min) / (vector.range[2] - vector.range[1])
  }

  (vector - vector.range[1]) * fac + min
}

data.paths <- function(analysis.path, actives = FALSE) {
    if (actives) {
        actives.path <- paste0(analysis.path, 'actives.csv')
    } else {
        actives.path <- NULL
    }

    matrix.data.path <- paste0(analysis.path, 'mat01.csv')
    selected.col <- paste0(analysis.path, 'selected.csv')

    list(matrix.path = matrix.data.path, selected = selected.col,
        actives.path = actives.path)
}

load.data <- function(parameters, analysis.path) {
    if (parameters$type != 'simimatrix' ||
        parameters$type == 'simiclustermatrix') {

        paths.list <- data.paths(analysis.path, actives = TRUE)

        matrix.data <- Matrix::readMM(paths.list$matrix.path)
        actives <- read.table(paths.list$actives.path, sep = '\t', quote = '"')

        colnames(matrix.data) <- actives[, 1]
    } else if (parameters$type == 'simimatrix' ||
        parameters$type == 'simiclustermatrix') {

        paths.list <- data.paths(analysis.path)

        matrix.data <- read.csv2(paths.list$matrix.path)
        matrix.data <- as.matrix(matrix.data)
    }

    list(paths = paths.list, matrix = matrix.data)
}

load.coords <- function() {
    coords <- try(coords, TRUE)
    load(paste0(analysis.path, 'RData.RData'))

    if (!is.matrix(coords)) {
        coords <- NULL
    }

    coords
}

select.matrix.word <- function(parameters, matrix.data, selected.column) {
    if ('word' %in% parameters) {
        word <- TRUE
        index <- parameters$word + 1
    } else {
        word <- FALSE
        index <- NULL
    }

    if (!word) {
        matrix.data <- matrix.data[, selected.column]
    } else {
        forme <- colnames(matrix.data)[index]

        if (!index %in% selected.column) {
            selected.column <- append(selected.column, index)
        }

        matrix.data <- matrix.data[, selected.column]
        index <- which(colnames(matrix.data) == forme)
    }

    list(matrix = matrix.data, index = index)
}

check.selected.column <- function(sel.col.file, matrix.data) {
    if (file.exists(sel.col.file)) {
        sel.col <- read.csv2(sel.col.file, header = FALSE)
        sel.col <- sel.col[, 1] + 1
    } else {
        sel.col <- 1:ncol(matrix.data)
    }

    sel.col
}

check.inf <- function(sparse) {
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

check.word.parameter <- function(parameters, sparse, matrix.data, index) {
    if ('word' %in% parameters) {
        sparse <- graph.word(sparse, index)
        col.sum <- colSums(sparse)

        if (length(col.sum)) sparse <- sparse[, -which(col.sum == 0)]

        row.sum <- rowSums(sparse)

        if (length(row.sum)) sparse <- sparse[-which(row.sum == 0),]
        if (length(col.sum)) matrix.data <- matrix.data[, -which(col.sum == 0)]

    }

    list(sparse = sparse, matrix = matrix.data)
}

create.sparse <- function(parameters, matrix.data) {
    if (parameters$method == 'cooc') {
        sparse <- square_matrix(matrix.data)

    } else if (parameters$method == 'Russel') {
        sparse <- proxy::simil(matrix.data, method = parameters$method,
            diag = TRUE, upper = TRUE, by_rows = FALSE)

    } else if (parameters$method == 'binomial') {
        sparse <- binom_sim(matrix.data)

    } else {
        sparse <- proxy::simil(as.matrix(matrix.data), method = parameters$method,
            diag = TRUE, upper = TRUE, by_rows = FALSE)
    }

    as.matrix(stats::as.dist(sparse, diag = TRUE, upper = TRUE))
}

check.coeff.tv <- function(parameters) {
    if (parameters$coeff.tv && parameters$sfromchi) {
        parameters$coeff.tv <- NULL
        parameters$minmaxeff <- c(parameters$tvmin, parameters$tvmax)
    } else {
        parameters$minmaxeff <- NULL
    }

    parameters
}

check.vcex <- function(parameters) {
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
#' @param analysis.path of data
#' @return list of objects
#' 
#' @export
generate.graph <- function(parameters, analysis.path) {
    parameters <- check.coeff.tv(parameters)
    parameters <- check.vcex(parameters)

    if (parameters$keep.coord) {
        parameters$coords <- load.coords()
    } else {
        parameters$coords <- NULL
    }

    data.definition <- load.data(parameters, analysis.path)

    selected.colmun <- check.selected.column(data.definition$paths$selected,
        data.definition$matrix)
    selection <- select.matrix.word(parameters, data.definition$matrix, selected.colmun)

    sparse <- create.sparse(parameters, selection$matrix)
    sparse <- check.inf(sparse)

    valid.data <- check.word.parameter(parameters, sparse, selection$matrix,
        selection$index)

    eff <- colSums(as.matrix(valid.data$matrix))
    x <- list(mat = valid.data$sparse, eff = eff)

    do_simi(
        x,
        method = parameters$method,
        seuil = parameters$seuil,
        p.type = parameters$type,
        layout.type = parameters$layout,
        max.tree = parameters$maxtree,
        coeff.vertex = parameters$coeff.tv,
        coeff.edge = parameters$coeff.te.range,
        minmaxeff = parameters$minmaxeff,
        vcexminmax = parameters$vcexminmax,
        cex = parameters$cex,
        coords = parameters$coords,
        communities = parameters$communities,
        halo = parameters$halo
    )
}

create.graph <- function(data) {
    data <- as.matrix(data)
    network <- igraph::graph_from_adjacency_matrix(data, mode = 'undirected')

    # Conta as arestas para cada par de vértices, atribui aos pesos das arestas
    # e realiza a simplificação do grafo (reduz a quantidade de aresta para 1)
    igraph::E(network)$weight <- igraph::count_multiple(network)
    network <- igraph::simplify(network, edge.attr.comb = 'max')
}
