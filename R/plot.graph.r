open_file_graph <- function (filename, width = 800, height = 800, quality = 100,
    svg = FALSE) {

	if (Sys.info()["sysname"] == 'Darwin') {
        width <- width / 74.97
        height <- height / 74.97

        if (!svg) {
		    quartz(file = filename, type = 'png', width = width, height = height)
        } else {
            svg(filename.to.svg(filename), width = width, height = height)
        }
	} else {
        if (svg) {
            svg(filename.to.svg(filename), width = width / 74.97,
                height = height / 74.97)
        } else {
		    png(filename, width = width, height = height)
        }
	}
}

apply.chd <- function(parameters, dm) {
    et <- list()

    for (index in seq_along(parameters$listet)) {
        line.et <- parameters$listet[[index]]
        line.et <- line.et + 1

        et[[index + 1]] <- paste(line.et, collapse = ',')
    }

    unetoile <- paste(parameters$selected.stars, collapse = "','")
    
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
    dcol <- apply(lex[[4]], 1, which.max)
    toblack <- apply(lex[[4]], 1, max)
    gcol <- rainbow(length(unetoile))
    
    vertex.label.color <- gcol[dcol]
    vertex.label.color[which(toblack <= 3.84)] <- 'black'
    
    leg <- list(unetoile = unetoile, gcol = gcol)
    parameters$cols <- vertex.label.color
    chi.vertex.size <- norm.vec(toblack, parameters$vcexmin,  parameters$vcexmax)
    
    list(to.black = toblack, vertex.label.color = vertex.label.color, leg = leg,
        chi.vertex.size = chi.vertex.size)
}

apply.plot.definitions <- function(parameters, graph.simi) {
    vertex.label.color <- 'black'
    chi.vertex.size <- 1
    leg <- NULL

    if ((parameters$type == 'clustersimitxt' &&
        parameters$tmpchi) ||
        parameters$type %in% c('simimatrix', 'simiclustermatrix') &&
        'tmpchi' %in% parameters) {
        
        lchi <- read.table(parameters$tmpchi)
        lchi <- lchi[, 1]
        lchi <- lchi[sel.col]
    }

    if (parameters$type %in% c('clustersimitxt', 'simimatrix', 'simiclustermatrix') &&
        parameters$cex.from.chi) {
        
        label.cex <- norm.vec(lchi, parameters$vcexmin, parameters$vcexmax)
    } else {
        if (is.null(parameters$vcexmin)) {
            label.cex <- parameters$cex
        } else {
            label.cex <- graph.simi$label.cex
        }
    }

    if (parameters$type %in% c('clustersimitxt', 'simimatrix', 'simiclustermatrix') &&
        parameters$sfromchi) {

        vertex.size <- norm.vec(lchi, parameters$tvmin, parameters$tvmax)
        
        if (!length(vertex.size)) vertex.size <- 0
    } else {
        if (is.null(parameters$tvmin)) {
            vertex.size <- 0
        } else {
            vertex.size <- graph.simi$eff
        }
    }

    list(vertex.size = vertex.size, leg = leg, chi.vertex.size = chi.vertex.size,
        vertex.label.color = vertex.label.color)
}

plot.graph <- function(graph.simi, parameters, ...) {
    if (parameters$bystar) {
        chd.definition <- apply.chd(parameters, matrix.data)
        vertex.label.color <- chd.definition$vertex.label.color
        leg <- chd.definition$leg

        if (parameters$cex.from.chi) {
            label.cex <- chd.definition$chi.vertex.size
        } else {
            label.cex <- parameters$cex
        }

        if (parameters$sfromchi) {
            vertex.size <- norm.vec(chd.definition$to.black, parameters$tvmin,
                parameters$tvmax)
        } else {
            vertex.size <- NULL
        }
    } else {
        plot.definitions <- apply.plot.definitions(parameters, graph.simi)
        vertex.size <- plot.definitions$vertex.size
        leg <- plot.definitions$leg
        label.cex <- plot.definitions$chi.vertex.size
        vertex.label.color <- plot.definitions$vertex.label.color
    }

    if (!is.null(graph.simi$communities)) {
        colm <- rainbow(length(graph.simi$communities))

        if (vertex.size != 0 || graph.simi$halo) {
            vertex.label.color <- 'black'
            parameters$cols <- colm[membership(graph.simi$communities)]
        } else {
            vertex.label.color <- colm[membership(graph.simi$communities)]
        }
    }

    coords <- plot.simi(
        graph.simi,
        ...,
        vertex.label = parameters$label.v,
        edge.label = parameters$label.e,
        vertex.col = parameters$cols,
        vertex.label.color = vertex.label.color,
        vertex.label.cex = label.cex,
        vertex.size = vertex.size,
        edge.col = parameters$cola,
        leg = leg,
        alpha = parameters$alpha,
        edge.curved = parameters$edge.curved,
        svg = parameters$svg
    )
}
