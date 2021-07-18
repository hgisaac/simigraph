test_that('generate.graph and plot.graph returns correct types', {
    analysis.path <- '../../iramuteq_results/corpus_simitxt/'

    parameters <- list(
        method = 'cooc', # 'Russel', 'binomial', etc.
        keep.coord = FALSE,
        type = 'simitxt',
        seuil = 0.01, # float, NULL
        layout = 'frutch', # 'random', 'circle', 'frutch', 'kawa', 'graphopt'
        maxtree = TRUE,
        coeff_vertex = 0,
        coeff.te.range = c(1, 10), # c(min, max), NULL
        coeff.tv = 0, # int, NULL
        sfromchi = FALSE,
        tvmin = 5,
        tvmax = 30,
        vcex = TRUE,
        vcexmin = 1.0,
        vcexmax = 2.5,
        cex.from.chi = FALSE,
        cex = 1.0,
        communities = NULL, # int, NULL
        halo = FALSE,
        
        # Used to plot
        alpha = 0.20,
        bystar = FALSE,
        svg = FALSE,
        cols = c(255, 0, 0, 255),
        cola = c(200, 200, 200, 255),
        label.v = 1,
        label.e = 1,
        edge.curved = TRUE
    )

    graph_simi <- generate.graph(parameters, analysis.path)

    expect_type(graph_simi, 'list')
    expect_type(plot.graph(graph_simi, parameters, p.type = 'nplot',
        filename = '../../images_test/graph_simi_refac.png'), 'double')
})
