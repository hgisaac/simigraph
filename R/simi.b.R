simiClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "simiClass",
    inherit = simiBase,
    private = list(
        .run = function() {
            if (length(self$data) == 0) return()

            result <- preprocess(self$data, 'corpus', min_docfreq = 4)
            
            graph_simi <- generate_graph(
                result$dtm,
                method = 'cooc',
                keep_coord = FALSE,
                seuil = 0.01,
                plot_type = 'nplot',
                layout_type = 'frutch',
                max_tree = TRUE,
                coeff_vertex = 0,
                coeff_edge_range = c(1, 10),
                sfromchi = FALSE,
                minmax_eff = c(5, 30),
                vcex_minmax = c(1.0, 2.5),
                cex = 1.0,
                communities = 1
            )
            
            self$results$plot$setState(graph_simi)
        },
        .plot = function(image, ...) {
            graph_simi <- image$state

            if (is.null(graph_simi)) return(FALSE)
            
            plot_result <- plot_graph(
                graph_simi,
                coeff_vertex = 0,
                cex_from_chi = FALSE,
                cex = 1.0,
                sfromchi = FALSE,
                minmax_eff = c(5, 30),
                vcex_minmax = c(1.0, 2.5),
                filename = 'graph_simi_refac.png',
                communities = graph_simi$communities,
                halo = TRUE
                #variable = 'codinome'
            )
            
            print(plot_result)
            TRUE
        }
    )
)
