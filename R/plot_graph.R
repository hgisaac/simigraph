open_file_graph <- function (
    filename,
    width = 800,
    height = 800,
    quality = 100,
    svg = FALSE
) {
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

apply_plot_definitions <- function(
    label_cex,
    eff,
    tmp_chi,
    cex_from_chi,
    vcex_minmax,
    cex,
    sfromchi,
    minmax_eff
) {
    vertex_label_color <- 'black'
    chi_vertex_size <- 1

    if (!is.null(tmp_chi)) {
        lchi <- read.table(tmp_chi)
        lchi <- lchi[, 1]
        lchi <- lchi[sel_col]
    }

    if (cex_from_chi) {
        label_cex <- norm_vec(lchi, vcex_minmax[1], vcex_minmax[2])
    } else {
        if (is.null(vcex_minmax)) {
            label_cex <- cex
        } else {
            label_cex <- label_cex
        }
    }

    if (sfromchi) {
        vertex_size <- norm_vec(lchi, minmax_eff[1], minmax_eff[2])
        
        if (!length(vertex_size)) vertex_size <- 0
    } else {
        if (is.null(minmax_eff)) {
            vertex_size <- 0
        } else {
            vertex_size <- eff
        }
    }

    list(
        vertex_size = vertex_size,
        chi_vertex_size = chi_vertex_size,
        vertex_label_color = vertex_label_color
    )
}

plot_graph <- function(
    graph_simi,
    coeff_vertex,
    cex_from_chi,
    cex,
    sfromchi,
    minmax_eff,
    vcex_minmax,
    halo = FALSE,
    plot_type = 'nplot',
    filename = NULL,
    communities = NULL,
    edge_color = 'black',
    vertex_label_color = 'black',
    vertex_label_cex = NULL,
    vertex_size = NULL,
    leg = NULL,
    width = 800,
    height = 800,
    alpha = 0.1,
    cexalpha = FALSE,
    edge_curved = TRUE,
    svg = FALSE,
    bg = 'white',
    tmp_chi = NULL,
    vertex_color = c(255, 0, 0, 255),
    variable = NULL
) {
    if (!is.null(variable)) {
        mapping <- map_variables(graph_simi$mat, variable, vcex_minmax)
        leg <- list(variables = mapping$var_values, colors = mapping$colors)

        if (cex_from_chi) {
            vertex_label_cex <- mapping$labels
        } else {
            vertex_label_cex <- cex
        }

        if (sfromchi) {
            vertex_size <- get_vertices_chi(graph_simi$mat, vcex_minmax)
        } else {
            vertex_size <- NULL
        }
    } else {
        plot_definitions <- apply_plot_definitions(
            vertex_label_cex,
            coeff_vertex,
            tmp_chi,
            cex_from_chi,
            vcex_minmax,
            cex,
            sfromchi,
            minmax_eff
        )
        vertex_size <- plot_definitions$vertex_size
        leg <- NULL
        vertex_label_cex <- plot_definitions$chi_vertex_size
        vertex_label_color <- plot_definitions$vertex_label_color
    }

    if (!is.null(communities)) {
        colm <- rainbow(length(communities))

        if (vertex_size != 0 || halo) {
            vertex_label_color <- 'black'
            vertex_color <- colm[igraph::membership(communities)]
        } else {
            vertex_label_color <- colm[igraph::membership(communities)]
        }
    }

    if (!is.null(vertex_label_cex)) {
        vertex_label_cex <- vertex_label_cex
    } else {
        vertex_label_cex <- graph_simi$label_cex
    }

    if (cexalpha) {
        vertex_label_color <- define_vertex_label_color(
            vertex_label_cex,
            vertex_label_color
        )
    }

    if (is.null(vertex_size)) {
        vertex_size <- graph_simi$coeff_vertex
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
            vertex_label_cex,
            edge_color,
            edge_curved,
            vertex_label_color,
            halo
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

map_variables <- function(dtm, variable, vcex_minmax) {
    doc_vars <- quanteda::docvars(dtm, variable)
    var_values <- unique(doc_vars)
    uces_variables <- list()
    
    for (index in seq_along(var_values)) {
        uces_variables[index] <- paste(
            which(doc_vars %in% var_values[index]),
            collapse = ','
        )
    }
    
    fsum <- NULL
    
    for (index in seq_along(var_values)) {
        tosum <- strsplit(uces_variables[[index]], ',')[[1]]

        if (length(tosum) > 1) {
            fsum <- cbind(fsum, Matrix::colSums(dtm[tosum, ]))
        } else {
            fsum <- cbind(fsum, dtm[tosum, ])
        }
    }

    labels <- map_labels(fsum, vcex_minmax)
    colors <- map_colors(fsum, var_values)

    list(labels = labels, colors = colors, var_values = var_values)
}

map_colors <- function(matrix_sum, variables) {
    # Create a vector of max chi values form matrix
    vertices_chi <- calculate_matrix_chi(matrix_sum)
    max_locations <- apply(vertices_chi, 1, which.max)

    # Query the color pallet from max chi values
    color_pallet <- RColorBrewer::brewer.pal(length(variables), 'Set1')
    vertices_color <- color_pallet[max_locations]

    # Assign lower chi vertices to 'non-variable'
    max_vertices <- apply(vertices_chi, 1, max)
    vertices_color[which(max_vertices <= 3.18)] <- 'black'
}

map_labels <- function(matrix_sum, vcex_minmax) {
    chi_matrix <- calculate_matrix_chi(matrix_sum)
    max_row <- apply(chi_matrix, 1, max)
    norm_vec(max_row, vcex_minmax[1],  vcex_minmax[2])
}

calculate_matrix_chi <- function(matrix_data) {
    matrix(
        sapply(seq_along(matrix_data), function(x) {
            result_matrix <- matrix(0, 2, 2)
            row_sum_matrix <- column_sum_matrix <- matrix_data
            
            row_sum_matrix[, 1:ncol(matrix_data)] <- rowSums(matrix_data)
            row_sum_matrix <- row_sum_matrix - matrix_data

            column_sum_matrix[1:nrow(matrix_data),] <- colSums(matrix_data)
            column_sum_matrix <- column_sum_matrix - matrix_data

            total_sum_matrix <- matrix(
                sum(rowSums(matrix_data)),
                nrow(matrix_data),
                ncol(matrix_data)
            )
            total_sum_matrix <- total_sum_matrix - row_sum_matrix - column_sum_matrix

            result_matrix[1, 1] <- matrix_data[x]
            result_matrix[1, 2] <- row_sum_matrix[x]
            result_matrix[2, 1] <- column_sum_matrix[x]
            result_matrix[2, 2] <- total_sum_matrix[x]

            chi_result <- chi_sq(result_matrix)
            
            if (is.na(chi_result$p_value)) {
                chi_result$p_value <- 1
                chi_result$statistic <- 0
            }

            if (result_matrix[1, 1] > chi_result$expected[1, 1]) {
                chi_result$statistic
            } else {
                0
            }
        }),
        ncol = ncol(matrix_data)
    )
}

chi_sq <- function(x) {
    expected <- outer(rowSums(x), colSums(x)) / sum(x)
    statistic <- sum(abs(x - expected) ^ 2 / expected)
    p_value <- stats::pchisq(statistic, 1, lower.tail = FALSE)
    list(statistic = statistic, expected = expected, p_value = p_value)
}

define_vertex_label_color <- function(vertex_label_cex, vertex_label_color) {
    alphas <- norm_vec(vertex_label_cex, 0.5, 1)
    new_vertex_label_color <- c()
    
    if (length(vertex_label_color) == 1) {
        for (i in seq_len(length(alphas))) {
            new_vertex_label_color <- append(
                new_vertex_label_color,
                adjustcolor(vertex_label_color, alpha = alphas[i])
            )
        }
    } else {
        for (i in seq_len(length(alphas))) {
            new_vertex_label_color <- append(
                new_vertex_label_color,
                adjustcolor(vertex_label_color[i], alpha = alphas[i])
            )
        }
    }

    new_vertex_label_color
}

n_plot <- function(
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
    vertex_label_color,
    halo
) {    
    open_file_graph(filename, width = width, height = height, svg = svg)
    par(mar = c(2, 2, 2, 2), bg = bg, pch = ' ')
    
    if (!is.null(leg)) {
        layout(
            matrix(c(1, 2), 1, 2, byrow = TRUE),
            widths = c(3, lcm(7))
        )
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
        if (halo) {
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
        legend(x = 'center', leg$variables, fill = leg$colors)
    }
    
    dev.off()
    graph_simi$layout
}

tk_plot <- function(
    graph_simi,
    vertex_size,
    vertex_color,
    vertex_label_color,
    edge_color
) { 
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
