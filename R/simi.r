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
    graph <- igraph::graph_from_adjacency_matrix(data_matrix, mode = 'lower',
        weighted = TRUE)
}

define_weights <- function(simi_graph, max.tree, method, seuil, mat.simi,
    minmaxeff, vcexminmax, coeff.vertex, cex, coeff.edge, mat.eff) {
    
    weori <- igraph::get.edge.attribute(simi_graph, 'weight')
    
    if (max.tree) {
        if (method == 'cooc') {
            invw <- 1 / weori
        } else {
            invw <- 1 - weori
        }

        igraph::E(simi_graph)$weight <- invw
        g.toplot <- igraph::minimum.spanning.tree(simi_graph)
        
        if (method == 'cooc') {
            igraph::E(g.toplot)$weight <- 1 / igraph::E(g.toplot)$weight
        } else {
            igraph::E(g.toplot)$weight <- 1 - igraph::E(g.toplot)$weight
        }
    }

    if (!is.null(seuil)) {
        if (seuil >= max(mat.simi)) seuil <- 0
        
        vec <- vector()
        tovire <- which(igraph::E(g.toplot)$weight <= seuil)
        g.toplot <- igraph::delete.edges(g.toplot, tovire)
        
        for (i in 1:(length(igraph::V(g.toplot)))) {
            if (length(igraph::neighbors(g.toplot, i)) == 0) {
                vec <- append(vec, i)
            }
        }

        g.toplot <- igraph::delete.vertices(g.toplot, vec)
        v.label <- igraph::V(g.toplot)$name

        if (!is.logical(vec)) mat.eff <- mat.eff[-vec]
    } else {
        v.label <- NULL
        vec <- NULL
    }

    if (!is.null(minmaxeff[1])) {
        eff <- norm.vec(mat.eff, minmaxeff[1], minmaxeff[2])
    } else {
        eff <- coeff.vertex
    }

    if (!is.null(vcexminmax[1])) {
        label.cex <- norm.vec(mat.eff, vcexminmax[1], vcexminmax[2])
    } else {
        label.cex <- cex
    }

    if (!is.null(coeff.edge)) {
        we.width <- norm.vec(abs(igraph::E(g.toplot)$weight), coeff.edge[1], coeff.edge[2])
    } else {
        we.width <- NULL
    }

    if (method != 'binom') {
        we.label <- round(igraph::E(g.toplot)$weight, 2)
    } else {
        we.label <- round(igraph::E(g.toplot)$weight, 3)
    }

    list(
        simi_graph = g.toplot,
        elim = vec,
        vertices_labels = v.label,
        eff = eff,
        label_cex = label.cex,
        we_width = we.width,
        we_label = we.label,
        mat_eff = mat.eff
    )
}

define_layout <- function(p.type, layout.type, coords, simi_graph) {
    layout_types <- list(
        frutch = igraph::layout.fruchterman.reingold,
        kawa = igraph::layout.kamada.kawai,
        random = igraph::layout.random,
        circle = igraph::layout.circle,
        graphopt = igraph::layout.graphopt
    )

    if (is.null(coords)) {
        layout_function <- layout_types[[layout.type]]

        if (!is.null(layout_function)) {
            defined_layout <- layout_function(simi_graph)
        }
    } else {
        defined_layout <- coords
    }

    defined_layout
}

define_communities <- function(communities, simi_graph) {
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

    communities_function <- communities_types[communities]
    
    if (!is.null(communities_function)) {
        communities_function(simi_graph)
    }
}

do_simi <- function(
    x,
    method = 'cooc',
    seuil = NULL,
    p.type = 'tkplot',
    layout.type = 'frutch',
    max.tree = TRUE,
    coeff.vertex = NULL,
    coeff.edge = NULL,
    minmaxeff = NULL,
    vcexminmax = NULL,
    cex = 1,
    coords = NULL,
    communities = NULL,
    halo = FALSE
) {
    simi_graph <- create_graph(x$mat)
    
    weights_definition <- define_weights(simi_graph, max.tree, method,
        seuil, x$mat, minmaxeff, vcexminmax, coeff.vertex, cex, coeff.edge, x$eff)

    graph_layout <- define_layout(p.type, layout.type, coords,
        weights_definition$simi_graph)

    if (!is.null(communities)) {
        graph_communities <- define_communities(communities,
            weights_definition$simi_graph)
    } else {
        graph_communities <- NULL
    }
    
    list(
        graph = weights_definition$simi_graph,
        mat.eff = weights_definition$mat_eff,
        eff = weights_definition$eff,
        mat = x$mat,
        halo = halo,
        layout = graph_layout,
        v.label = weights_definition$vertices_labels,
        we.width = weights_definition$we_width,
        we.label = weights_definition$we_label,
        communities = graph_communities,
        label.cex = weights_definition$label_cex,
        elim = weights_definition$elim
    )
}

define_vertex_label_color <- function(vertex_label_cex, vertex_label_color) {
    alphas <- norm.vec(vertex_label_cex, 0.5, 1)
    new_vertex_label_color <- c()
    
    if (length(vertex_label_color) == 1) {
        for (i in 1:length(alphas)) {
            new_vertex_label_color <- append(new_vertex_label_color,
                adjustcolor(vertex_label_color, alpha = alphas[i]))
        }
    } else {
        for (i in 1:length(alphas)) {
            new_vertex_label_color <- append(new_vertex_label_color,
                adjustcolor(vertex_label_color[i], alpha = alphas[i]))
        }
    }

    new_vertex_label_color
}

n_plot <- function(filename, width, height, svg, bg, leg, graph_simi, vertex.size,
    vertex.color, label.cex, edge.color, edge.curved, vertex_label_color) {
    
    open_file_graph(filename, width = width, height = height, svg = svg)
    par(mar = c(2, 2, 2, 2))
    par(bg = bg)
    
    if (!is.null(leg)) {
        layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(3, lcm(7)))
        par(mar = c(2, 2, 1, 0))
    }

    par(pch = ' ')

    if (is.null(graph_simi$com)) {
        plot(
            graph_simi$graph,
            vertex.label = '',
            edge.width = graph_simi$we.width,
            vertex.size = vertex.size,
            vertex.color = vertex.color,
            vertex.label.color = 'white',
            edge.label = graph_simi$we.label,
            edge.label.cex = label.cex,
            edge.color = edge.color,
            vertex.label.cex = 0,
            layout = graph_simi$layout,
            edge.curved = edge.curved
        )
    } else {
        if (graph_simi$halo) {
            mark.groups <- communities(graph_simi$com)
        } else {
            mark.groups <- NULL
        }

        plot(
            graph_simi$com,
            graph_simi$graph,
            vertex.label = '',
            edge.width = graph_simi$we.width,
            vertex.size = vertex.size,
            vertex.color = vertex.color,
            vertex.label.color = 'white',
            edge.label = graph_simi$we.label,
            edge.label.cex = label.cex,
            edge.color = edge.color,
            vertex.label.cex = 0,
            layout = graph_simi$layout,
            mark.groups = mark.groups,
            edge.curved = edge.curved
        )
    }
    
    txt.layout <- igraph::layout.norm(graph_simi$layout, -1, 1, -1, 1, -1, 1)
    text(
        txt.layout[, 1],
        txt.layout[, 2],
        graph_simi$v.label,
        cex = label.cex,
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

tk_plot <- function(graph_simi, vertex.size, vertex.color, vertex_label_color,
    edge.color) {
    
    id <- tkplot(
        graph_simi$graph,
        vertex.label = graph_simi$v.label,
        edge.width = graph_simi$we.width,
        vertex.size = vertex.size,
        vertex.color = vertex.color,
        vertex.label.color = vertex_label_color,
        edge.label = graph_simi$we.label,
        edge.color = edge.color,
        layout = graph_simi$layout
    )

    coords <- tkplot.getcoords(id)
    ok <- try(
        coords <- tkplot.getcoords(id),
        TRUE
    )
    
    while (is.matrix(ok)) {
        ok <- try(
            coords <- tkplot.getcoords(id),
            TRUE
        )
        Sys.sleep(0.5)
    }
    
    tkplot.off()
    coords
}

plot_simi <- function(
    graph.simi,
    p.type = 'tkplot',
    filename = NULL,
    communities = NULL,
    vertex.color = 'red',
    edge.color = 'black',
    vertex.label.color = 'black',
    vertex.label.cex = NULL,
    vertex.size = NULL,
    leg = NULL,
    width = 800,
    height = 800,
    alpha = 0.1,
    cexalpha = FALSE,
    movie = NULL,
    edge.curved = TRUE,
    svg = FALSE,
    bg = 'white'
) {
    if (!is.null(vertex.label.cex)) {
        label.cex <- vertex.label.cex
    } else {
        label.cex <- graph.simi$label.cex
    }

    if (cexalpha) {
        vertex.label.color <- define_vertex_label_color(label.cex, vertex.label.color)
    }

    if (is.null(vertex.size)) {
        vertex.size <- graph.simi$eff
    }

    if (p.type == 'nplot') {
        plot_result <- n_plot(
            filename,
            width,
            height,
            svg,
            bg,
            leg,
            graph.simi,
            vertex.size,
            vertex.color,
            label.cex,
            edge.color,
            edge.curved,
            vertex.label.color
        )
    } else if (p.type == 'tkplot') {
        plot_result <- tk_plot(
            graph.simi,
            vertex.size,
            vertex.color,
            vertext.label.color,
            edge.color
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
        res$p.value <- 0
    }

    res$p.value
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
    
    for (i in 1:(ncol(x) - 1)) {
        for (j in (i + 1):ncol(x)) {
            mat[j, i] <- make_bin(cs, a, i, j, n)
        }
    }

    mat
}

graph.word <- function(mat.simi, index) {
    nm <- matrix(
        0,
        ncol = ncol(mat.simi),
        nrow = nrow(mat.simi),
        dimnames = list(row.names(mat.simi), colnames(mat.simi))
    )

    nm[, index] <- mat.simi[, index]
    nm[index,] <- mat.simi[index,]
}
