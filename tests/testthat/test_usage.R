parameters <- list(
    method = 'cooc', # 'Russel', 'binomial', etc.
    keep_coord = FALSE,
    type = 'simitxt',
    seuil = 0.01, # float, NULL
    plot_type = 'nplot',
    layout = 'frutch', # 'random', 'circle', 'frutch', 'kawa', 'graphopt'
    max_tree = TRUE,
    coeff_vertex = 0,
    coeff_edge_range = c(1, 10), # c(min, max), NULL
    coeff_tv = 0, # int, NULL
    sfromchi = FALSE,
    tvmin = 5,
    tvmax = 30,
    vcex = TRUE,
    vcex_min = 1.0,
    vcex_max = 2.5,
    cex_from_chi = FALSE,
    cex = 1.0,
    communities = NULL, # int, NULL
    halo = FALSE,
    
    # Used to plot
    alpha = 0.20,
    bystar = FALSE,
    svg = FALSE,
    cols = c(255, 0, 0, 255),
    cola = c(200, 200, 200, 255),
    label_v = 1,
    label_e = 1,
    edge_curved = TRUE
)

result <- preprocess('../../data_test/corpus.csv', 'corpus', min_docfreq = 4)
graph_simi <- generate_graph(parameters, result$dtm)

test_that('generate_graph returns list', {
    expect_type(graph_simi, 'list')
})

test_that('plot_graph returns double', {
    expect_type(
        plot_graph(
            graph_simi,
            parameters,
            filename = '../../images_test/graph_simi_refac.png'
        ),
        'double'
    )
})
