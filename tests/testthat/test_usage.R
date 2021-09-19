data <- read.csv('../../data_test/corpus.csv')
result <- preprocess(data, 'corpus', min_docfreq = 4)
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

test_that('generate_graph returns list', {
    expect_type(graph_simi, 'list')
})

test_that('plot_graph returns double', {
    expect_type(
        plot_graph(
            graph_simi,
            coeff_vertex = 0,
            cex_from_chi = FALSE,
            cex = 1.0,
            sfromchi = FALSE,
            minmax_eff = c(5, 30),
            vcex_minmax = c(1.0, 2.5),
            filename = '../../images_test/graph_simi_refac.png',
            communities = graph_simi$communities,
            halo = TRUE,
            variable = 'codinome'
        ),
        'double'
    )
})
