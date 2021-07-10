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

define_weights <- function(simi_graph, max.tree, method, weori, seuil, mat.simi,
    minmaxeff, vcexminmax, coeff.vertex, cex, coeff.edge, mat.eff) {
    
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

    list(simi_graph = g.toplot, elim = vec, vertices_labels = v.label, eff = eff,
        label_cex = label.cex, we_width = we.width, we_label = we.label, mat_eff = mat.eff)
}

define_layout <- function(p.type, layout.type, coords, simi_graph) {
    if (p.type == 'rgl' || p.type == 'rglweb') {
        nd <- 3
    } else {
        nd <- 2
    }

    if (is.null(coords)) {
        if (layout.type == 'frutch')
            lo <- igraph::layout.fruchterman.reingold(simi_graph, dim = nd)
        
        if (layout.type == 'kawa')
            lo <- igraph::layout.kamada.kawai(simi_graph, dim = nd)
        
        if (layout.type == 'random')
            lo <- igraph::layout.random(simi_graph, dim = nd)
        
        if (layout.type == 'circle' & p.type != 'rgl')
            lo <- igraph::layout.circle(simi_graph)
        
        if (layout.type == 'circle' & p.type == 'rgl')
            lo <- igraph::layout.sphere(simi_graph)
        
        if (layout.type == 'graphopt')
            lo <- igraph::layout.graphopt(simi_graph)
    } else {
        lo <- coords
    }

    lo
}

define_communities <- function(communities, simi_graph) {
    if (communities == 0 ){
        graph_communities <- igraph::edge.betweenness.community(simi_graph)
    } else if (communities == 1) {
        graph_communities <- igraph::fastgreedy.community(simi_graph)
    } else if (communities == 2) {
        graph_communities <- igraph::label.propagation.community(simi_graph)
    } else if (communities == 3) {
        graph_communities <- igraph::leading.eigenvector.community(simi_graph)
    } else if (communities == 4) {
        graph_communities <- igraph::multilevel.community(simi_graph)
    } else if (communities == 5) {
        graph_communities <- igraph::optimal.community(simi_graph)
    } else if (communities == 6) {
        graph_communities <- igraph::spinglass.community(simi_graph)
    } else if (communities == 7) {
        graph_communities <- igraph::walktrap.community(simi_graph)
    }

    graph_communities
}

do.simi <- function(x, method = 'cooc', seuil = NULL, p.type = 'tkplot',
    layout.type = 'frutch', max.tree = TRUE, coeff.vertex = NULL, coeff.edge = NULL,
    minmaxeff = NULL, vcexminmax = NULL, cex = 1, coords = NULL, communities = NULL,
    halo = FALSE) {
    
    simi_graph <- create_graph(x$mat)
    weori <- igraph::get.edge.attribute(simi_graph, 'weight')
    
    weights_definition <- define_weights(simi_graph, max.tree, method, weori,
        seuil, x$mat, minmaxeff, vcexminmax, coeff.vertex, cex, coeff.edge, x$eff)

    graph_layout <- define_layout(p.type, layout.type, coords,
        weights_definition$simi_graph)

    if (!is.null(communities)) {
        graph_communities <- define_communities(communities,
            weights_definition$simi_graph)
    } else {
        graph_communities <- NULL
    }
    
    list(graph = weights_definition$simi_graph, mat.eff = weights_definition$mat_eff,
        eff = weights_definition$eff, mat = x$mat, halo = halo, layout = graph_layout,
        v.label = weights_definition$vertices_labels, we.width = weights_definition$we_width,
        we.label = weights_definition$we_label, communities = graph_communities,
        label.cex = weights_definition$label_cex, elim = weights_definition$elim)
}

plot.simi <- function(graph.simi, p.type = 'tkplot',filename = NULL, communities = NULL,
    vertex.col = 'red', edge.col = 'black', edge.label = TRUE, vertex.label = TRUE,
    vertex.label.color = 'black', vertex.label.cex= NULL, vertex.size=NULL, leg = NULL,
    width = 800, height = 800, alpha = 0.1, cexalpha = FALSE, movie = NULL,
    edge.curved = TRUE, svg = FALSE, bg = 'white') {
    
    mat.simi <- graph.simi$mat
    g.toplot <- graph.simi$graph
    if (is.null(vertex.size)) {
        vertex.size <- graph.simi$eff
    } else {
        vertex.size <- vertex.size
    }
    we.width <- graph.simi$we.width
    if (vertex.label) {
        #v.label <- vire.nonascii(graph.simi$v.label)
        v.label <- graph.simi$v.label
    } else {
        v.label <- NA
    }
    if (edge.label) {
        we.label <- graph.simi$we.label
    } else {
        we.label <- NA
    }
    lo <- graph.simi$layout
    if (!is.null(vertex.label.cex)) {
        label.cex<-vertex.label.cex
    } else {
        label.cex = graph.simi$label.cex
    }
    if (cexalpha) {
        alphas <- norm.vec(label.cex, 0.5,1)
        nvlc <- NULL
        if (length(vertex.label.color) == 1) {
            for (i in 1:length(alphas)) {
             nvlc <- append(nvlc, adjustcolor(vertex.label.color, alpha=alphas[i]))
            }
        } else {
            for (i in 1:length(alphas)) {
                nvlc <- append(nvlc, adjustcolor(vertex.label.color[i], alpha=alphas[i]))
            }
        }
        vertex.label.color <- nvlc  
    }
    if (p.type=='nplot') {
        #print('ATTENTION - PAS OPEN FILE')
        open_file_graph(filename, width = width, height = height, svg = svg)
        par(mar=c(2,2,2,2))
        par(bg=bg)
        
        if (!is.null(leg)) {
            layout(matrix(c(1,2),1,2, byrow=TRUE),widths=c(3,lcm(7)))
            par(mar=c(2,2,1,0))
        }

        par(pch=' ')

        if (is.null(graph.simi$com)) {
            plot(g.toplot, vertex.label = '', edge.width = we.width, vertex.size = vertex.size,
                vertex.color = vertex.col, vertex.label.color = 'white', edge.label = we.label,
                edge.label.cex = label.cex, edge.color = edge.col, vertex.label.cex = 0,
                layout = lo, edge.curved = edge.curved)
        } else {
            if (graph.simi$halo) {
                mark.groups <- communities(graph.simi$com)
            } else {
                mark.groups <- NULL
            }
            plot(graph.simi$com, g.toplot,vertex.label = '', edge.width = we.width,
                vertex.size = vertex.size, vertex.color=vertex.col, vertex.label.color='white',
                edge.label=we.label, edge.label.cex=label.cex, edge.color=edge.col,
                vertex.label.cex = 0, layout=lo, mark.groups = mark.groups, edge.curved=edge.curved)
        }
        #txt.layout <- lo
        txt.layout <- igraph::layout.norm(lo, -1, 1, -1, 1, -1, 1)
        #txt.layout <- txt.layout[order(label.cex),]
        #vertex.label.color <- vertex.label.color[order(label.cex)]
        #v.label <- v.label[order(label.cex)]
        #label.cex <- label.cex[order(label.cex)]
        text(txt.layout[,1], txt.layout[,2], v.label, cex=label.cex, col=vertex.label.color)
        if (!is.null(leg)) {
            par(mar=c(0,0,0,0))
            plot(0, axes = FALSE, pch = '')
            legend(x = 'center' , leg$unetoile, fill = leg$gcol)
        }
        dev.off()
        return(lo)
    }
    if (p.type=='tkplot') {
        id <- tkplot(g.toplot,vertex.label=v.label, edge.width=we.width, vertex.size=vertex.size,
            vertex.color=vertex.col, vertex.label.color=vertex.label.color, edge.label=we.label,
            edge.color=edge.col, layout=lo)

        coords <- tkplot.getcoords(id)
        ok <- try(coords <- tkplot.getcoords(id), TRUE)
        while (is.matrix(ok)) {
            ok <- try(coords <- tkplot.getcoords(id), TRUE)
            Sys.sleep(0.5)
        }
    tkplot.off()
    return(coords)
    }
    
    if (p.type == 'rgl' || p.type == 'rglweb') {
        library('rgl')
        #rgl.open()
        #par3d(cex=0.8)
        lo <- layout.norm(lo, -10, 10, -10, 10, -10, 10)
        bg3d('white')
        rglplot(g.toplot,vertex.label='', edge.width=we.width/10, vertex.size=0.01,
            vertex.color=vertex.col, vertex.label.color="black", edge.color = edge.col,
            layout=lo, rescale = FALSE)
        
        text3d(lo[,1], lo[,2], lo[,3], vire.nonascii(v.label), col = vertex.label.color,
            alpha = 1, cex = vertex.label.cex)

        rgl.spheres(lo, col = vertex.col, radius = vertex.size/100, alpha = alpha)
        #rgl.bg(color = c('white','black'))
        #bg3d('white')
        if (!is.null(movie)) {
            require(tcltk)
            ReturnVal <- tkmessageBox(title="RGL 3 D", message = "Cliquez pour commencer le film",
                icon="info",type="ok")

            movie3d(spin3d(axis=c(0,1,0),rpm=6), movie = 'film_graph', frames = "tmpfilm",
                duration=10, clean=TRUE, top = TRUE, dir = movie)

            ReturnVal <- tkmessageBox(title="RGL 3 D",message="Film fini !",
                icon="info",type="ok")
        }
        
        if (p.type == 'rglweb') {
            writeWebGL(dir = filename, width = width, height= height)
        } else {
            require(tcltk)
            ReturnVal <- tkmessageBox(title="RGL 3 D",message="Cliquez pour fermer",
                icon="info",type="ok")
        }

        rgl.close()
    }
}

make.bin <- function(cs, a, i, j, nb) {
    if (a[i, j] >= 1) {
        ab <- a[i, j] - 1 
        res <- binom.test(ab, nb, (cs[i]/nb) * (cs[j]/nb), "less")
    } else {
        res <- NULL
        res$p.value <- 0
    }

    res$p.value
}

binom.sim <- function(x) {
    a <- t(x) %*% (x)
    n <- nrow(x)
    cs <- colSums(x)
    mat <- matrix(0,ncol(x),ncol(x))
    colnames(mat)<-colnames(a)
    rownames(mat)<-rownames(a)
    for (i in 1:(ncol(x)-1)) {
        for (j in (i+1):ncol(x)) {
            mat[j,i] <- make.bin(cs, a, i, j , n)
        }
    }

    mat
}

graph.word <- function(mat.simi, index) {
    nm <- matrix(0, ncol = ncol(mat.simi), nrow=nrow(mat.simi),
        dimnames = list(row.names(mat.simi), colnames(mat.simi)))

    nm[,index] <- mat.simi[,index]
    nm[index,] <- mat.simi[index,]
    nm
}
