test_that('preprocess returns uce_uc_table and dtm', {
    result <- preprocess('../../data_test/corpus.csv', 'corpus')
    expect_s3_class(result$uce_uc_table, 'data.frame')
    expect_s4_class(result$dtm, 'dfm')
})
