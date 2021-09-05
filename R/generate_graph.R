create_graph <- function(data_matrix, mode = 'lower', weighted = TRUE, ...) {
    igraph::graph_from_adjacency_matrix(
        data_matrix,
        mode = mode,
        weighted = weighted,
        ...
    )
}

invert_weights <- function(simi_graph, method) {
    if (method == 'cooc') {
        igraph::E(simi_graph)$weight <- 1 / igraph::E(simi_graph)$weight
    } else {
        igraph::E(simi_graph)$weight <- 1 - igraph::E(simi_graph)$weight
    }
    
    simi_graph
}

define_weights <- function(simi_graph, method) {
    simi_graph <- invert_weights(simi_graph, method)
    simi_graph <- igraph::minimum.spanning.tree(simi_graph)
    invert_weights(simi_graph, method)
}

simplify_graph <- function(simi_graph, seuil, mat_simi, mat_eff) {
    if (!is.null(seuil)) {
        if (seuil >= max(mat_simi)) seuil <- 0
        
        vec <- vector()
        tovire <- which(igraph::E(simi_graph)$weight <= seuil)
        simi_graph <- igraph::delete.edges(simi_graph, tovire)
        
        for (i in seq_len(length(igraph::V(simi_graph)))) {
            if (length(igraph::neighbors(simi_graph, i)) == 0) {
                vec <- append(vec, i)
            }
        }

        simi_graph <- igraph::delete.vertices(simi_graph, vec)

        if (!is.logical(vec)) mat_eff <- mat_eff[-vec]
    } else {
        vec <- NULL
    }

    list(
        simi_graph = simi_graph,
        elim = vec,
        mat_eff = mat_eff
    )
}

get_labels <- function(simi_graph, seuil, method) {
    if (!is.null(seuil)) {
        vertice_label <- igraph::V(simi_graph)$name
    } else {
        vertice_label <- NULL
    }

    if (method != 'binom') {
        digits <- 2
    } else {
        digits <- 3
    }

    edge_label <- round(igraph::E(simi_graph)$weight, digits)

    list(vertices = vertice_label, edges = edge_label)
}

define_width <- function(simi_graph, coeff_edge) {
    if (!is.null(coeff_edge)) {
        edge_width <- norm_vec(
            abs(igraph::E(simi_graph)$weight),
            coeff_edge[1],
            coeff_edge[2]
        )
    } else {
        edge_width <- NULL
    }

    edge_width
}

define_layout <- function(plot_type, layout_type, coords, simi_graph) {
    layout_types <- list(
        frutch = igraph::layout.fruchterman.reingold,
        kawa = igraph::layout.kamada.kawai,
        random = igraph::layout.random,
        circle = igraph::layout.circle,
        graphopt = igraph::layout.graphopt
    )

    if (is.null(coords)) {
        layout_function <- layout_types[[layout_type]]

        if (!is.null(layout_function)) {
            defined_layout <- layout_function(simi_graph)
        }
    } else {
        defined_layout <- coords
    }

    defined_layout
}

define_communities <- function(communities_number, simi_graph) {
    communities_types <- c(
        igraph::edge.betweenness.community,
        igraph::fastgreedy.community,
        igraph::label.propagation.community,
        igraph::leading.eigenvector.community,
        igraph::multilevel.community,
        igraph::optimal.community,
        igraph::spinglass.community,
        igraph::walktrap.community
    )

    communities_function <- communities_types[[communities_number]]
    
    if (!is.null(communities_function)) {
        return(communities_function(simi_graph))
    }
}

make_bin <- function(cs, a, i, j, nb) {
    if (a[i, j] >= 1) {
        ab <- a[i, j] - 1 
        res <- binom.test(ab, nb, (cs[i] / nb) * (cs[j] / nb), "less")
    } else {
        res <- NULL
        res$p_value <- 0
    }

    res$p_value
}

square_matrix <- function(x) {
    x <- as.matrix(x)
    t(x) %*% x
}

binom_sim <- function(x) {
    a <- square_matrix(x)
    n <- nrow(x)
    cs <- colSums(x)
    mat <- matrix(0, ncol(x), ncol(x))
    colnames(mat) <- colnames(a)
    rownames(mat) <- rownames(a)
    
    for (i in seq_along(ncol(x) - 1)) {
        for (j in (i + 1):ncol(x)) {
            mat[j, i] <- make_bin(cs, a, i, j, n)
        }
    }

    mat
}

graph_word <- function(mat_simi, index) {
    nm <- matrix(
        0,
        ncol = ncol(mat_simi),
        nrow = nrow(mat_simi),
        dimnames = list(row.names(mat_simi), colnames(mat_simi))
    )

    nm[, index] <- mat_simi[, index]
    nm[index,] <- mat_simi[index,]
}

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

check_word_parameter <- function(word, sparse, matrix_data) {
    sparse <- graph_word(sparse, word)
    col_sum <- colSums(sparse)

    if (length(col_sum)) sparse <- sparse[, -which(col_sum == 0)]

    row_sum <- rowSums(sparse)

    if (length(row_sum)) sparse <- sparse[-which(row_sum == 0),]
    if (length(col_sum)) matrix_data <- matrix_data[, -which(col_sum == 0)]

    list(sparse = sparse, matrix = matrix_data)
}

compute_similarity <- function(method, matrix_data) {
    if (method == 'cooc') {
        sparse <- square_matrix(matrix_data)

    } else if (method == 'Russel') {
        sparse <- proxy::simil(
            matrix_data,
            method = method,
            diag = TRUE,
            upper = TRUE,
            by_rows = FALSE
        )

    } else if (method == 'binomial') {
        sparse <- binom_sim(matrix_data)

    } else {
        sparse <- proxy::simil(
            as.matrix(matrix_data),
            method = method,
            diag = TRUE,
            upper = TRUE,
            by_rows = FALSE
        )
    }

    as.matrix(stats::as.dist(sparse, diag = TRUE, upper = TRUE))
}

#' Generate similitude graph
#' 
#' @param parameters list
#' @param analysis_path of data
#' 
#' @return list of objects
#' @export
generate_graph <- function(
    dtm,
    method = 'cooc',
    seuil = NULL,
    plot_type = 'nplot',
    layout_type = 'frutch',
    max_tree = TRUE,
    coeff_vertex = NULL,
    coeff_edge_range = NULL,
    minmax_eff = NULL,
    vcex_minmax = NULL,
    cex = 1.0,
    coords = NULL,
    communities = NULL,
    keep_coord = FALSE,
    sfromchi = FALSE,
    cex_from_chi = FALSE,
    word = NULL
) {
    if (coeff_vertex && sfromchi) {
        coeff_vertex <- NULL
    } else {
        minmax_eff <- NULL
    }
    
    if (cex_from_chi && is.null(vcex_minmax)) {
        stop('A range of value to vcex_minmax is necessary.')
    }

    if (keep_coord) {
        coords <- load_coords()
    }

    sparse <- compute_similarity(method, dtm)
    sparse <- check_inf(sparse)

    if (!is.null(word)) {
        valid_data <- check_word_parameter(word, sparse, dtm)
        sparse <- valid_data$sparse
        dtm <- valid_data$matrix
    }

    eff <- colSums(as.matrix(dtm))
    simi_graph <- create_graph(sparse)

    if (max_tree) {
        simi_graph <- define_weights(simi_graph, method)
    }
    
    simplified_graph <- simplify_graph(
        simi_graph,
        seuil,
        sparse,
        eff
    )
    
    graph_labels <- get_labels(
        simplified_graph$simi_graph,
        seuil,
        method
    )
    
    edge_width <- define_width(
        simplified_graph$simi_graph,
        coeff_edge_range
    )

    if (!is.null(minmax_eff)) {
        coeff_vertex <- norm_vec(
            eff,
            minmax_eff[1],
            minmax_eff[2]
        )
    }

    if (!is.null(vcex_minmax)) {
        cex <- norm_vec(
            eff,
            vcex_minmax[1],
            vcex_minmax[2]
        )
    }

    graph_layout <- define_layout(
        plot_type,
        layout_type,
        coords,
        simplified_graph$simi_graph
    )

    if (!is.null(communities)) {
        communities <- define_communities(
            communities,
            simplified_graph$simi_graph
        )
    }
    
    list(
        graph = simplified_graph$simi_graph,
        mat_eff = simplified_graph$mat_eff,
        coeff_vertex = coeff_vertex,
        mat = dtm,
        layout = graph_layout,
        v_label = graph_labels$vertices,
        we_width = edge_width,
        we_label = graph_labels$edges,
        communities = communities,
        label_cex = cex,
        elim = simplified_graph$elim
    )
}
