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
