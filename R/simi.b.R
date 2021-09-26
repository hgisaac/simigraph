simiClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "hcaClass",
    inherit = hcaBase,
    private = list(
        .run = function() {
            # code
        }
    )
)
