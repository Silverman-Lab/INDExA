#' @export
setGeneric("getReads", function(.object) standardGeneric("getReads"))
#' @export
setGeneric("getDirichletSamples", function(.object) standardGeneric("getDirichletSamples"))
#' @export
setGeneric("getLog2NormSamples", function(.object) standardGeneric("getLog2NormSamples"))
#' @export
setGeneric("getDenom", function(.object) standardGeneric("getDenom"))
#' @export
setGeneric("getScales", function(.object) standardGeneric("getScales"))
#' @rdname indexa.mc.function
#' @export
setGeneric("indexa.mc", function(reads, mc.samples=250, denom="none", custom.scale=NULL) standardGeneric("indexa.mc"), signature=c("reads"))
