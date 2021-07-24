open_file_graph <- function (filename, width = 800, height = 800, quality = 100,
    svg = FALSE) {

	if (Sys.info()["sysname"] == 'Darwin') {
        width <- width / 74.97
        height <- height / 74.97

        if (!svg) {
		    quartz(file = filename, type = 'png', width = width, height = height)
        } else {
            svg(filename_to_svg(filename), width = width, height = height)
        }
	} else {
        if (svg) {
            svg(filename_to_svg(filename), width = width / 74.97,
                height = height / 74.97)
        } else {
		    png(filename, width = width, height = height)
        }
	}
}

apply_chd <- function(parameters, dm) {
    et <- list()

    for (index in seq_along(parameters$listet)) {
        line_et <- parameters$listet[[index]]
        line_et <- line_et + 1

        et[[index + 1]] <- paste(line_et, collapse = ',')
    }

    unetoile <- paste(parameters$selected_stars, collapse = "','")
    
    fsum <- NULL
    rs <- rowSums(dm)
    
    for (i in 1:length(unetoile)) {
        print(unetoile[i])
        tosum <- et[[i]]

        if (length(tosum) > 1) {
            fsum <- cbind(fsum, colSums(dm[tosum, ]))
        } else {
            fsum <- cbind(fsum, dm[tosum, ])
        }
    }

    source('~/iramuteq-0.7-alpha2/Rscripts/chdfunct.r')

    lex <- AsLexico2(fsum, chip = TRUE)
    dcol <- apply(lex[[4]], 1, which_max)
    toblack <- apply(lex[[4]], 1, max)
    gcol <- rainbow(length(unetoile))
    
    vertex_label_color <- gcol[dcol]
    vertex_label_color[which(toblack <= 3.84)] <- 'black'
    
    leg <- list(unetoile = unetoile, gcol = gcol)
    parameters$cols <- vertex_label_color
    chi_vertex_size <- norm_vec(toblack, parameters$vcexmin,  parameters$vcexmax)
    
    list(to_black = toblack, vertex_label_color = vertex_label_color, leg = leg,
        chi_vertex_size = chi_vertex_size)
}

apply_plot_definitions <- function(parameters, graph_simi) {
    vertex_label_color <- 'black'
    chi_vertex_size <- 1
    leg <- NULL

    if ((parameters$type == 'clustersimitxt' &&
        parameters$tmpchi) ||
        parameters$type %in% c('simimatrix', 'simiclustermatrix') &&
        'tmpchi' %in% parameters) {
        
        lchi <- read.table(parameters$tmpchi)
        lchi <- lchi[, 1]
        lchi <- lchi[sel_col]
    }

    if (parameters$type %in% c('clustersimitxt', 'simimatrix', 'simiclustermatrix') &&
        parameters$cex_from_chi) {
        
        label_cex <- norm_vec(lchi, parameters$vcexmin, parameters$vcexmax)
    } else {
        if (is.null(parameters$vcexmin)) {
            label_cex <- parameters$cex
        } else {
            label_cex <- graph_simi$label_cex
        }
    }

    if (parameters$type %in% c('clustersimitxt', 'simimatrix', 'simiclustermatrix') &&
        parameters$sfromchi) {

        vertex_size <- norm_vec(lchi, parameters$tvmin, parameters$tvmax)
        
        if (!length(vertex_size)) vertex_size <- 0
    } else {
        if (is.null(parameters$tvmin)) {
            vertex_size <- 0
        } else {
            vertex_size <- graph_simi$eff
        }
    }

    list(vertex_size = vertex_size, leg = leg, chi_vertex_size = chi_vertex_size,
        vertex_label_color = vertex_label_color)
}

plot_graph <- function(graph_simi, parameters, ...) {
    if (parameters$bystar) {
        chd_definition <- apply_chd(parameters, matrix_data)
        vertex_label_color <- chd_definition$vertex_label_color
        leg <- chd_definition$leg

        if (parameters$cex_from_chi) {
            label_cex <- chd_definition$chi_vertex_size
        } else {
            label_cex <- parameters$cex
        }

        if (parameters$sfromchi) {
            vertex_size <- norm_vec(chd_definition$to_black, parameters$tvmin,
                parameters$tvmax)
        } else {
            vertex_size <- NULL
        }
    } else {
        plot_definitions <- apply_plot_definitions(parameters, graph_simi)
        vertex_size <- plot_definitions$vertex_size
        leg <- plot_definitions$leg
        label_cex <- plot_definitions$chi_vertex_size
        vertex_label_color <- plot_definitions$vertex_label_color
    }

    if (!is.null(graph_simi$communities)) {
        colm <- rainbow(length(graph_simi$communities))

        if (vertex_size != 0 || graph_simi$halo) {
            vertex_label_color <- 'black'
            parameters$cols <- colm[igraph::membership(graph_simi$communities)]
        } else {
            vertex_label_color <- colm[igraph::membership(graph_simi$communities)]
        }
    }

    coords <- plot_simi(
        graph_simi,
        ...,
        plot_type = parameters$plot_type,
        vertex_color = parameters$cols,
        vertex_label_color = vertex_label_color,
        vertex_label_cex = label_cex,
        vertex_size = vertex_size,
        edge_color = parameters$cola,
        leg = leg,
        alpha = parameters$alpha,
        edge_curved = parameters$edge_curved,
        svg = parameters$svg
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
