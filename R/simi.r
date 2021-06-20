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

define_weights <- function(graph) {
    if (max.tree) {
        if (method == 'cooc') {
		    invw <- 1 / weori
        } else {
            invw <- 1 - weori
        }
		E(g1)$weight<-invw
		g.max<-minimum.spanning.tree(g1)
        if (method == 'cooc') {
		    E(g.max)$weight<-1 / E(g.max)$weight
        } else {
            E(g.max)$weight<-1 - E(g.max)$weight
        }
		g.toplot<-g.max
	}

    if (!is.null(seuil)) {
        if (seuil >= max(mat.simi)) seuil <- 0
        vec<-vector()
        w<-E(g.toplot)$weight
        tovire <- which(w<=seuil)
        g.toplot <- delete.edges(g.toplot,(tovire))
        for (i in 1:(length(V(g.toplot)))) {
            if (length(neighbors(g.toplot,i))==0) {
                vec<-append(vec,i)
            }
        }
        g.toplot <- delete.vertices(g.toplot,vec)
        v.label <- V(g.toplot)$name
        if (!is.logical(vec)) mat.eff <- mat.eff[-(vec)]
    } else {
		vec <- NULL
	}

	if (!is.null(minmaxeff[1])) {
        eff<-norm.vec(mat.eff,minmaxeff[1],minmaxeff[2])
    } else {
        eff<-coeff.vertex
    }
    if (!is.null(vcexminmax[1])) {
        label.cex = norm.vec(mat.eff, vcexminmax[1], vcexminmax[2])
    } else {
        label.cex = cex
    }
    if (!is.null(coeff.edge)) {
        we.width <- norm.vec(abs(E(g.toplot)$weight), coeff.edge[1], coeff.edge[2]) 
	    #we.width <- abs((E(g.toplot)$weight/max(abs(E(g.toplot)$weight)))*coeff.edge)
    } else {
        we.width <- NULL
    }
    if (method != 'binom') {
        we.label <- round(E(g.toplot)$weight,2)
    } else {
        we.label <- round(E(g.toplot)$weight,3)
    }
}

define_layout <- function() {
    if (p.type=='rgl' || p.type=='rglweb') {
        nd<-3
    } else {
        nd<-2
    }

    if (is.null(coords)) {
    	if (layout.type == 'frutch')
    		lo <- layout.fruchterman.reingold(g.toplot,dim=nd)
    	if (layout.type == 'kawa')
    		lo <- layout.kamada.kawai(g.toplot,dim=nd)
    	if (layout.type == 'random')
    		lo <- layout.random(g.toplot,dim=nd)
    	if (layout.type == 'circle' & p.type != 'rgl')
    		lo <- layout.circle(g.toplot)
    	if (layout.type == 'circle' & p.type == 'rgl')
    		lo <- layout.sphere(g.toplot)
        if (layout.type == 'graphopt')
            lo <- layout.graphopt(g.toplot)
    } else {
        lo <- coords
    }
}

define_communities <- function() {
    if (communities == 0 ){ #'edge.betweenness.community') {
        com <- edge.betweenness.community(g.toplot)
    } else if (communities == 1) {
        com <- fastgreedy.community(g.toplot)
    } else if (communities == 2) {
        com <- label.propagation.community(g.toplot)
    } else if (communities == 3) {
        com <- leading.eigenvector.community(g.toplot)
    } else if (communities == 4) {
        com <- multilevel.community(g.toplot)
    } else if (communities == 5) {
        com <- optimal.community(g.toplot)
    } else if (communities == 6) {
        com <- spinglass.community(g.toplot)
    } else if (communities == 7) {
        com <- walktrap.community(g.toplot)
    }
}

do.simi <- function(x, method = 'cooc', seuil = NULL, p.type = 'tkplot',
    layout.type = 'frutch', max.tree = TRUE, coeff.vertex = NULL, coeff.edge = NULL,
    minmaxeff = NULL, vcexminmax = NULL, cex = 1, coords = NULL, communities = NULL,
    halo = FALSE) {
	
    mat.simi <- x$mat
    mat.eff <- x$eff
    v.label <- colnames(mat.simi)
	g1 <- create_graph(mat.simi)
	g.toplot <- g1
	weori <- get.edge.attribute(g1, 'weight')
	
    define_weights()
    
	define_layout()

    if (communities) {
        define_communities()
    } else {
        com <- NULL
    }
    
	list(graph = g.toplot, mat.eff = mat.eff, eff = eff, mat = mat.simi,
        v.label = v.label, we.width = we.width, we.label = we.label, label.cex = label.cex,
        layout = lo, communities = com, halo = halo, elim = vec)
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
        txt.layout <- layout.norm(lo, -1, 1, -1, 1, -1, 1)
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
