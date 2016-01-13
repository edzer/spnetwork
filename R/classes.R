#' @import methods sp igraph
#' @importFrom stats weights

#  @exportClass igraph
setClass("igraph") # S4

#' The SpatialNetwork class and constructor function
#'
#' SpatialNetwork: a Class for Spatial Networks (lines and edges)
#'
#' A class to store spatial networks, such as a road network.
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{bbox}}{matrix, holding bounding box}
#'    \item{\code{proj4string}}{object of class \link[sp]{CRS-class}}
#'    \item{\code{lines}}{list}
#'    \item{\code{data}}{data.frame, holding attributes associated with lines}
#'    \item{\code{g}}{object of a subclass of \link[igraph]{igraph}}
#'    \item{\code{nb}}{neighbourhood list}
#'    \item{\code{weightfield}}{character; describing the weight field used}
#'  }
#'
#' @usage SpatialNetwork(sl, g, nb, weights, weightfield)
#' @param sl object of one of (a sublasses of) \link{SpatialLines}
#' @param g object of class \link{igraph} 
#' @param nb neighbourhood list
#' @param weights weight for edges (defaults to length of linestring)
#' @param weightfield character; name of the attribute field of \code{sl} that will be used as weights
#' @return object of class \link{SpatialNetwork-class}

#' @name SpatialNetwork-class
#' @rdname SpatialNetwork-class
#' @aliases SpatialNetwork SpatialNetwork-class [[,SpatialNetwork,ANY,missing-method [[<-,SpatialNetwork,ANY,missing-method $,SpatialNetwork-method
#' @exportClass SpatialNetwork
#' @author Edzer Pebesma
#' @note note
#'
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
#' sl = SpatialLines(l)
#' sln = SpatialNetwork(sl)
#' plot(sln@g$x, sln@g$y, col = sln@g$n, pch = 16, cex = 2, asp = 1)
#' lines(sl)
#' library(igraph) # E
#' text(sln@g$x, sln@g$y, E(sln@g), pos = 4)
#' plot(sln@g)
setClass("SpatialNetwork",
	contains = "SpatialLinesDataFrame", 
	slots = c(g = "igraph", nb = "list", weightfield = "character"),
    validity = function(object) {
		# print("Entering validation: SpatialNetwork")
		if (any(sapply(object@lines, function(x) length(x@Lines)) != 1))
			stop("all Lines objects need to have a single Line") 
    	if (length(object@lines) != length(E(object@g)))
			stop("edges do not match number line segments") 
    	if (length(object@nb) != length(V(object@g)))
			stop("vertices do not match length of neighbourhood list") 
		return(TRUE)
	}
)

#' @export
SpatialNetwork = function(sl, g, nb, weights, weightfield) {
    stopifnot(is(sl, "SpatialLines"))
    if (!is(sl, "SpatialLinesDataFrame")) 
        sl = new("SpatialLinesDataFrame", sl, data = data.frame(id = 1:length(sl)))
    if (!all(sapply(sl@lines, function(x) length(x@Lines)) == 1)) 
        stop("SpatialLines is not simple: each Lines element should have only a single Line")
	if (missing(g) || missing(nb)) { # sort out from sl
    	startEndPoints = function(x) {
        	firstLast = function(L) {
            	cc = coordinates(L)[[1]]
            	rbind(cc[1, ], cc[nrow(cc), ])
        	}
        	do.call(rbind, lapply(x@lines, firstLast))
    	}
    	s = startEndPoints(sl)
    	zd = zerodist(SpatialPoints(s))
    	pts = 1:nrow(s)
	
    	# the following can't be done vector-wise, there is a progressive effect:
    	if (nrow(zd) > 0) {
        	for (i in 1:nrow(zd)) 
				pts[zd[i, 2]] = pts[zd[i, 1]]
    	}
    	stopifnot(identical(s, s[pts, ]))
	
    	# map to 1:length(unique(pts))
    	pts0 = match(pts, unique(pts))
    	node = rep(1:length(sl), each = 2)
    	g = graph(pts0, directed = FALSE)  # edges
    	# nb = lapply(1:length(unique(pts)), function(x) node[which(pts0 == x)])
		nb = lapply(as.list(incident_edges(g, V(g))), as.numeric) # meaning, we can drop it as slot
    	nodes = s[unique(pts), ]
    	g$x = nodes[, 1]  # x-coordinate vertex
    	g$y = nodes[, 2]  # y-coordinate vertex
    	g$n = as.vector(table(pts0))  # nr of edges
    	# line lengths:
    	sl$length = SpatialLinesLengths(sl) # takes care of projected
    	if (!missing(weights)) {
    	    if (missing(weightfield)) {
    	        weightfield = 'weight'
    	    }
	        sl@data[, weightfield] <- weights
    	} else
    	    weightfield = 'length'
        E(g)$weight = sl@data[, weightfield]
        
    	# create list with vertices, starting/stopping for each edge?  add for
    	# each SpatialLines, the start and stop vertex
    	pts2 = matrix(pts0, ncol = 2, byrow = TRUE)
    	sl$start = pts2[, 1]
    	sl$end = pts2[, 2]
	}
	if (!missing(weights)) {
	    if (missing(weightfield))
	        weightfield = 'weight'
        sl@data[, weightfield] <- weights
    	E(g)$weight = weights
	}
    new("SpatialNetwork", sl, g = g, nb = nb, weightfield = weightfield)
}

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


#' extract method for testClass
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
	sl = as(x, "SpatialLinesDataFrame")
	sl = sl[i, j, ..., drop = FALSE]
	V(x@g)$sel = 1:length(V(x@g)) # identify
	g = subgraph.edges(x@g, i)
	sel = V(g)$sel
	g$x = g$x[sel]
	g$y = g$y[sel]
    g$n = g$n[sel]
	nb = lapply(as.list(incident_edges(g, V(g))), as.numeric)
	new("SpatialNetwork", sl, g = g, nb = nb, weightfield = x@weightfield)
})

#' @param col color
#' @param cex symbol size
#' @rdname SpatialNetwork-methods
#' @export
points.SpatialNetwork = function(x, ..., col = "red", cex = 2) {
	points(x@g$x, x@g$y, col = col, cex = cex, ...)
}
