#' @export 
encodeSequences <- function(sequences) 
# Encodes base sequences into a 2-bit encoding.
{
    out <- .Call(cxx_encode_sequences, sequences)
    names(out) <- names(sequences)
    out
}
