#' Get or set weight field in SpatialNetwork
#'
#' Get or set value of weight field in SpatialNetwork
#'
#' @section Details:
#' These functions manipulate the value of weightfield in a
#' SpatialNetwork. When changing the value of weightfield, the weights
#' of the graph network are updated with the values of the corresponding
#' variables.
#'
#' @usage weights(object, ...)
#' @usage weigths(x) <- value
#' @param object SpatialNetwork to use
#' @param ... ignored
#' @param x SpatialNetwork to use
#' @param weightfield The name of the variable to set/use.
#' @param value Either the name of the variable to use as the weight field or
#' a vector containing the weights to use if \code{varname} is
#' passed to the replacement function.
#' @name weightfield
#' @aliases weights weights<- weights,SpatialNetwork-method weights<-,SpatialNetwork,character-method, weights<-,SpatialNetwork,vector-method, weights<-,SpatialNetwork,character,vector-method
#' @examples
#' library(sp) # Lines, SpatialLines
#' l0 = cbind(c(1, 2), c(0, 0))
#' l1 = cbind(c(0, 0, 0), c(0, 1, 2))
#' l2 = cbind(c(0, 0, 0), c(2, 3, 4))
#' l3 = cbind(c(0, 1, 2), c(2, 2, 3))
#' l4 = cbind(c(0, 1, 2), c(4, 4, 3))
#' l5 = cbind(c(2, 2), c(0, 3))
#' l6 = cbind(c(2, 3), c(3, 4))
#' LL = function(l, ID) Lines(list(Line(l)), ID)
#' l = list(LL(l0, "e"), LL(l1, "a"), LL(l2, "b"), LL(l3, "c"), LL(l4, "d"), LL(l5, "f"), LL(l6, "g"))
#' sln = SpatialNetwork(SpatialLines(l))
#' weights(sln)
#' weights(sln) = 2 * sln$length
#' weights(sln) = "length"
#' weights(sln, "randomweights") = runif(nrow(sln))
NULL

##' @export
#setGeneric("weights",
#           function(x) standardGeneric("weights"))

#' @export
setGeneric("weights<-",
           function(x, value) standardGeneric("weights<-"))

#' @export
setGeneric("weights<-",
           function(x, weightfield, value) standardGeneric("weights<-"))

#' @rdname weightfield
#' @export weights
#' @export
weights.SpatialNetwork = function(object, ...) {
    message(paste0("Weightfield = ", object@weightfield))
    object@data[, object@weightfield]
}

#' @rdname weightfield
setReplaceMethod("weights", signature(x = "SpatialNetwork", value = "character"), 
                 definition = function(x, value) {
                     if (value %in% names(x)) {
                         x@weightfield <- value
                         E(x@g)$weight <- x@data[,value]
                         x
                     } else
                         stop(paste0("No field of name ",value," - weights not updated"))
                 })

#' @rdname weightfield
setReplaceMethod("weights", signature(x = "SpatialNetwork", value = "vector"),
                 definition = function(x, value) {
                     message("Using 'weight' as field name.")
                     x@data[,'weight'] <- value
                     x@weightfield <- 'weight'
                     E(x@g)$weight <- x@data[,'weight']
                     x
                 })

#' @rdname weightfield
setReplaceMethod("weights", signature(x = "SpatialNetwork", weightfield = "character", value = "vector"),
                 definition = function(x, weightfield, value) {
                     x@data[,weightfield] <- value
                     x@weightfield <- weightfield
                     E(x@g)$weight <- x@data[,weightfield]
                     x
                 })


#' extract method for SpatialNetwork
#'
#' @name [
#' @aliases [,SpatialNetwork-method
#' @docType methods
#' @param x object of class SpatialNetwork
#' @param i numeric; features to select
#' @param j numeric or character; attributes to select
#' @param ... ignored
#' @param drop logical; ignored
#' @rdname SpatialNetwork-methods
setMethod("[", c("SpatialNetwork", "ANY", "ANY"), function(x, i, j, ... , drop = TRUE) {
	# select i: edge_ids
	sl = as(x, "SpatialLinesDataFrame")[i, j, ..., drop = FALSE]
	g = subgraph.edges(x@g, i)
	new("SpatialNetwork", sl, g = g, weightfield = x@weightfield)
})

#' @param col color
#' @param cex symbol size
#' @rdname SpatialNetwork-methods
#' @export
points.SpatialNetwork = function(x, ..., col = "red", cex = 2) {
	points(V(x@g)$x, V(x@g)$y, col = col, cex = cex, ...)
}

if (!isGeneric("summary"))
	setGeneric("summary", function(object, ...)
		standardGeneric("summary"))

#' summary method for SpatialNetwork
#'
#' @name summary
#' @aliases summary,SpatialNetwork-method
#' @docType methods
#' @param object object of class SpatialNetwork
#' @rdname SpatialNetwork-methods
summary.SpatialNetwork = function(object, ...) {
    obj = list()
	obj[["class"]] = class(object)
    obj[["bbox"]] = bbox(object)
    obj[["is.projected"]] = is.projected(object)
    obj[["proj4string"]] = object@proj4string@projargs
	if ("data" %in% slotNames(object) && ncol(object@data) > 0)
		obj[["data"]] = summary(object@data)
	obj$edges = length(object)
	obj$nodes = length(V(object@g))
	obj$weightfield = object@weightfield
	obj$g = object@g
    class(obj) = "summary.SpatialNetwork"
    obj
}
setMethod("summary", "SpatialNetwork", summary.SpatialNetwork)

#' @rdname SpatialNetwork-methods
#' @export
print.summary.SpatialNetwork = function(x, ...) {
    cat(paste("Object of class ", x[["class"]], "\n", sep = ""))
    cat("Coordinates:\n")
    print(x[["bbox"]], ...)
    cat(paste("Is projected:", x[["is.projected"]], "\n"))
#    cat(paste("proj4string : [", x[["proj4string"]], "]\n", sep=""))
    pst <- paste(strwrap(x[["proj4string"]]), collapse="\n")
    if (nchar(pst) < 40) cat(paste("proj4string : [", pst, "]\n", sep=""))
    else cat(paste("proj4string :\n[", pst, "]\n", sep=""))
	cat(paste("# edges:", x$edges, "# nodes/vertices:", x$nodes, "\n"))
	cat(paste("# weightfield:", x$weightfield, "\n"))
	cat("Graph summary:\n")
	summary(x$g) # directly prints
    if (!is.null(x$data)) {
        cat("Lines attributes:\n")
        print(x$data, ...)
    }
    invisible(x)
}
