data <- as.matrix(read.csv('../../data_test/coocurrences.csv', row.name = 1))
simi_graph <- create_graph(data)

test_that('create_graph returns igraph object', {
    expect_s3_class(simi_graph, 'igraph')
})

test_that('define_communities returns communities object', {
    expect_s3_class(define_communities(1, simi_graph), 'communities')
})
