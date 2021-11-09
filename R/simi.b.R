simiClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "simiClass",
    inherit = simiBase,
    private = list(
        .run = function() {
            if (length(self$data) == 0) return()

            data <- self$data
            data[[self$options$text]] <- as.character(
                data[[self$options$text]]
            )

            result <- preprocess(
                data,
                self$options$text,
                self$options$min_seg_size,
                segment_size = self$options$seg_size,
                language = self$options$language,
                min_docfreq = self$options$min_docfreq
            )
            
            graph_simi <- generate_graph(
                result$dtm,
                method = self$options$method,
                seuil = self$options$seuil,
                plot_type = 'nplot',
                layout_type = self$options$layout_type,
                coeff_vertex = self$options$coeff_vertex,
                coeff_edge_range = c(
                    self$options$min_coeff_edge,
                    self$options$max_coeff_edge
                ),
                sfromchi = self$options$size_from_chi,
                minmax_eff = c(
                    self$options$min_eff,
                    self$options$max_eff
                ),
                vcex_minmax = c(
                    self$options$min_vcex,
                    self$options$max_vcex
                ),
                cex = self$options$cex,
                communities = self$options$communities
            )
            
            self$results$plot$setState(graph_simi)
        },
        .plot = function(image, ...) {
            graph_simi <- image$state

            if (is.null(graph_simi)) return(FALSE)
            
            plot_graph(
                graph_simi,
                coeff_vertex = self$options$coeff_vertex,
                cex_from_chi = self$options$cex_from_chi,
                cex = self$options$cex,
                sfromchi = self$options$size_from_chi,
                minmax_eff = c(
                    self$options$min_eff,
                    self$options$max_eff
                ),
                vcex_minmax = c(
                    self$options$min_vcex,
                    self$options$max_vcex
                ),
                communities = graph_simi$communities,
                halo = self$options$halo,
                variable = self$options$variable
            )

            TRUE
        }
    )
)
