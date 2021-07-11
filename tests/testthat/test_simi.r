test_that('create_graph returns igraph object', {
    data <- as.matrix(read.csv('../../data_test/coocurrences.csv', row.name = 1))
    expect_s3_class(create_graph(data), 'igraph')
})

test_that('do.simi is a closure', {
    expect_type(do_simi, 'closure')
})
