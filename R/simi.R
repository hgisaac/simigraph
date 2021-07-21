# a, b, c, and d are the counts of all (TRUE, TRUE), (TRUE, FALSE), (FALSE, TRUE), and (FALSE, FALSE)
# n <- a + b + c + d = nrow(x)

# jaccard a, b, c   a / (a + b + c)
# Kulczynski1 a, b, c   a / (b + c)
# Kulczynski2 a, b, c   [a / (a + b) + a / (a + c)] / 2
# Mountford a, b, c    2a / (ab + ac + 2bc)
# Fager, McGowan a, b, c   a / sqrt((a + b)(a + c)) - 1 / 2 sqrt(a + c)
# Russel, Rao a (a/n)
# Dice, Czekanowski, Sorensen a, b, c   2a / (2a + b + c)
# Mozley, Margalef a, b, c  an / (a + b)(a + c)
# Ochiai a, b, c  a / sqrt[(a + b)(a + c)]
# Simpson a, b, c   a / min{(a + b), (a + c)}
# Braun-Blanquet a, b, c  a / max{(a + b), (a + c)}

# Simple matching
# Sokal/Michener a, b, c, d, ((a + d) /n)
# Hamman, a, b, c, d, ([a + d] - [b + c]) / n
# Faith , a, b, c, d, (a + d/2) / n
# Tanimoto, Rogers a, b, c, d, (a + d) / (a + 2b + 2c + d)
# Phi  a, b, c, d   (ad - bc) / sqrt[(a + b)(c + d)(a + c)(b + d)]
# Stiles a, b, c, d  log(n(|ad-bc| - 0.5n)^2 / [(a + b)(c + d)(a + c)(b + d)])
# Michael   a, b, c, d   4(ad - bc) / [(a + d)^2 + (b + c)^2]
# Yule a, b, c, d  (ad - bc) / (ad + bc)
# Yule2  a, b, c, d  (sqrt(ad) - sqrt(bc)) / (sqrt(ad) + sqrt(bc))

create_graph <- function(data_matrix) {
    igraph::graph_from_adjacency_matrix(data_matrix, mode = 'lower',
        weighted = TRUE)
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

do_simi <- function(
    x,
    method = 'cooc',
    seuil = NULL,
    plot_type = 'nplot',
    layout_type = 'frutch',
    max_tree = TRUE,
    coeff_vertex = NULL,
    coeff_edge = NULL,
    minmax_eff = NULL,
    vcex_minmax = NULL,
    cex = 1,
    coords = NULL,
    communities = NULL,
    halo = FALSE
) {
    simi_graph <- create_graph(x$mat)

    if (max_tree) {
        simi_graph <- define_weights(simi_graph, method)
    }
    
    simplified_graph <- simplify_graph(simi_graph, seuil,
        x$mat, x$eff)
    
    graph_labels <- get_labels(simplified_graph$simi_graph, seuil, method)
    
    edge_width <- define_width(simplified_graph$simi_graph, coeff_edge)

    if (!is.null(minmax_eff[1])) {
        coeff.vertex <- norm_vec(x$eff, minmax_eff[1], minmax_eff[2])
    }

    if (!is.null(vcex_minmax[1])) {
        cex <- norm_vec(x$eff, vcex_minmax[1], vcex_minmax[2])
    }

    graph_layout <- define_layout(plot_type, layout_type, coords,
        simplified_graph$simi_graph)

    if (!is.null(communities)) {
        communities <- define_communities(communities,
            simplified_graph$simi_graph)
    }
    
    list(
        graph = simplified_graph$simi_graph,
        mat_eff = simplified_graph$mat_eff,
        eff = coeff_vertex,
        mat = x$mat,
        halo = halo,
        layout = graph_layout,
        v_label = graph_labels$vertices,
        we_width = edge_width,
        we_label = graph_labels$edges,
        communities = communities,
        label_cex = cex,
        elim = simplified_graph$elim
    )
}

define_vertex_label_color <- function(vertex_label_cex, vertex_label_color) {
    alphas <- norm_vec(vertex_label_cex, 0.5, 1)
    new_vertex_label_color <- c()
    
    if (length(vertex_label_color) == 1) {
        for (i in seq_len(length(alphas))) {
            new_vertex_label_color <- append(new_vertex_label_color,
                adjustcolor(vertex_label_color, alpha = alphas[i]))
        }
    } else {
        for (i in seq_len(length(alphas))) {
            new_vertex_label_color <- append(new_vertex_label_color,
                adjustcolor(vertex_label_color[i], alpha = alphas[i]))
        }
    }

    new_vertex_label_color
}

n_plot <- function(filename, width, height, svg, bg, leg, graph_simi, vertex_size,
    vertex_color, label_cex, edge_color, edge_curved, vertex_label_color) {
    
    open_file_graph(filename, width = width, height = height, svg = svg)
    par(mar = c(2, 2, 2, 2), bg = bg, pch = ' ')
    
    if (!is.null(leg)) {
        layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(3, lcm(7)))
        par(mar = c(2, 2, 1, 0))
    }

    if (is.null(graph_simi$communities)) {
        plot(
            graph_simi$graph,
            vertex.label = '',
            edge.width = graph_simi$we_width,
            vertex.size = vertex_size,
            vertex.color = vertex_color,
            vertex.label.color = 'white',
            edge.label = graph_simi$we_label,
            edge.label.cex = label_cex,
            edge.color = edge_color,
            vertex.label.cex = 0,
            layout = graph_simi$layout,
            edge.curved = edge_curved
        )
    } else {
        if (graph_simi$halo) {
            mark_groups <- igraph::communities(graph_simi$communities)
        } else {
            mark_groups <- NULL
        }

        plot(
            graph_simi$communities,
            graph_simi$graph,
            vertex.label = '',
            edge.width = graph_simi$we_width,
            vertex.size = vertex_size,
            vertex.color = vertex_color,
            vertex.label.color = 'white',
            edge.label = graph_simi$we_label,
            edge.label.cex = label_cex,
            edge.color = edge_color,
            vertex.label.cex = 0,
            layout = graph_simi$layout,
            mark.groups = mark_groups,
            edge.curved = edge_curved
        )
    }
    
    txt_layout <- igraph::norm_coords(graph_simi$layout)
    text(
        txt_layout[, 1],
        txt_layout[, 2],
        graph_simi$v_label,
        cex = label_cex,
        col = vertex_label_color
    )
    
    if (!is.null(leg)) {
        par(mar = c(0, 0, 0, 0))
        plot(0, axes = FALSE, pch = '')
        legend(x = 'center', leg$unetoile, fill = leg$gcol)
    }
    
    dev.off()
    graph_simi$layout
}

tk_plot <- function(graph_simi, vertex_size, vertex_color, vertex_label_color,
    edge_color) {
    
    id <- igraph::tkplot(
        graph_simi$graph,
        vertex.label = graph_simi$v_label,
        edge.width = graph_simi$we_width,
        vertex.size = vertex_size,
        vertex.color = vertex_color,
        vertex.label.color = vertex_label_color,
        edge.label = graph_simi$we_label,
        edge.color = edge_color,
        layout = graph_simi$layout
    )

    coords <- igraph::tkplot.getcoords(id)
    ok <- try(
        coords <- igraph::tkplot.getcoords(id),
        TRUE
    )
    
    while (is.matrix(ok)) {
        ok <- try(
            coords <- igraph::tkplot.getcoords(id),
            TRUE
        )
        Sys.sleep(0.5)
    }
    
    igraph::tkplot.off()
    coords
}

plot_simi <- function(
    graph_simi,
    plot_type = 'nplot',
    filename = NULL,
    communities = NULL,
    vertex_color = 'red',
    edge_color = 'black',
    vertex_label_color = 'black',
    vertex_label_cex = NULL,
    vertex_size = NULL,
    leg = NULL,
    width = 800,
    height = 800,
    alpha = 0.1,
    cexalpha = FALSE,
    movie = NULL,
    edge_curved = TRUE,
    svg = FALSE,
    bg = 'white'
) {
    if (!is.null(vertex_label_cex)) {
        label_cex <- vertex_label_cex
    } else {
        label_cex <- graph_simi$label_cex
    }

    if (cexalpha) {
        vertex_label_color <- define_vertex_label_color(label_cex, vertex_label_color)
    }

    if (is.null(vertex_size)) {
        vertex_size <- graph_simi$eff
    }

    if (plot_type == 'nplot') {
        plot_result <- n_plot(
            filename,
            width,
            height,
            svg,
            bg,
            leg,
            graph_simi,
            vertex_size,
            vertex_color,
            label_cex,
            edge_color,
            edge_curved,
            vertex_label_color
        )
    } else if (plot_type == 'tkplot') {
        plot_result <- tk_plot(
            graph_simi,
            vertex_size,
            vertex_color,
            vertex_label_color,
            edge_color
        )
    }
    plot_result
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
