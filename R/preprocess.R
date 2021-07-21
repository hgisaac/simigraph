preprocess <- function(
    file_path,
    column,
    min_segment_size = 0,
    doc_id = NULL,
    segment_size = 40,
    language = 'en',
    ...
) {
    corpus <- create_corpus(file_path, column)
    corpus <- rainette::split_segments(corpus, segment_size = segment_size)
    
    dtm <- clean_corpus(corpus, language, ...)
    dtm <- rainette::merge_segments(dtm, min_segment_size = min_segment_size,
        doc_id = doc_id)
    
    uce_uc_table <- create_uce_uc_table(dtm)
    dtm <- group_by_uc(dtm)

    list(uce_uc_table = uce_uc_table, dtm = dtm)
}

create_corpus <- function(file_path, column) {
    file_data <- read.csv(file_path)
    quanteda::corpus(file_data, text_field = column)
}

clean_corpus <- function(corpus, language, ...) {
    tokens <- quanteda::tokens(corpus, remove_punct = TRUE)
    tokens <- quanteda::tokens_remove(tokens, stopwords::stopwords(language))
    
    dtm <- quanteda::dfm(tokens, tolower = TRUE)
    dtm <- quanteda::dfm_trim(dtm, ...)
    dtm <- quanteda::dfm_remove(dtm, '')
    quanteda::dfm_weight(dtm, scheme = 'boolean')
}

create_uce_uc_table <- function(dtm) {
    data.frame(
        uce = quanteda::docvars(dtm, 'rainette_uce_id'),
        uc = quanteda::docvars(dtm, 'rainette_uc_id')
    )
}

group_by_uc <- function(dtm) {
    dtm <- quanteda::dfm_group(dtm, quanteda::docvars(dtm, 'rainette_uc_id'))
    quanteda::dfm_weight(dtm, scheme = 'boolean', force = TRUE)
}
