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
setReplaceMethod("weights", signature(x = "SpatialNetwork", weightfield = "missing", 
	value = "character"), 
                 definition = function(x, value) {
                     if (value %in% names(x)) {
                         x@weightfield <- value
                         E(x@g)$weight <- x@data[,value]
                         x
                     } else
                         stop(paste0("No field of name ",value," - weights not updated"))
                 })

#' @rdname weightfield
setReplaceMethod("weights", signature(x = "SpatialNetwork", weightfield = "missing", 
	value = "vector"),
                 definition = function(x, value) {
                     message("Using 'weight' as field name.")
                     x@data[,'weight'] <- value
                     x@weightfield <- 'weight'
                     E(x@g)$weight <- x@data[,'weight']
                     x
                 })

#' @rdname weightfield
setReplaceMethod("weights", signature(x = "SpatialNetwork", weightfield = "character", 
	value = "vector"),
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
#' @param i numeric; Lines to select (and any edges referring to them)
#' @param j numeric; graph edges to select (and any Lines they refer to)
#' @param ... ignored
#' @param drop logical; ignored
#' @rdname SpatialNetwork-methods
setMethod("[", c("SpatialNetwork", "ANY", "ANY"), function(x, i, j, ... , drop = TRUE) {
	if (!missing(i) && !missing(j))
		stop("use obj[i,] for selecting lines, obj[,j] for selecing edges, but not both")
	# select i: Lines ids
	if (!missing(i)) {
		if (length(i) > length(x) || any(i < 1 | i > length(x)))
			stop("index i out of range (1 ... # lines)")
		sl = as(x, "SpatialLinesDataFrame")[i, j, ..., drop = FALSE]
		if (is_directed(x@g))
			i = which(E(x@g)$link_index %in% i)
		g = subgraph.edges(x@g, i)
		E(g)$link_index = match(i, unique(i)) # map to 1..length(sl)
	} else if (!missing(j)) {
		g = subgraph.edges(x@g, j)
		sl = as(x, "SpatialLinesDataFrame")[unique(E(g)$link_index), ]
	} else
		return(x)
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

#' summary method for SpatialNetwork
#'
#' @name summary
#' @aliases summary,SpatialNetwork-method
#' @param object object of class SpatialNetwork
#' @docType methods
#' @rdname SpatialNetwork-methods
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

#' spTransform method for SpatialNetwork
#'
#' @name spTransform
#' @aliases spTransform,SpatialNetwork,ANY-method
#' @docType methods
#' @param CRSobj object of class \link[sp]{CRS}
#' @rdname SpatialNetwork-methods
setMethod("spTransform", signature("SpatialNetwork", "ANY"),
    function(x, CRSobj, ...) {
        new("SpatialNetwork",
			spTransform(as(x, "SpatialLinesDataFrame"), CRSobj, ...), # calls the rgdal methods
			g = x@g, weightfield = x@weightfield)
    }
)

#' As("SpatialNetwork", "SpatialPointsDataFrame")
#'
#' @aliases coerce,SpatialNetwork,SpatialPointsDataFrame-method
#' @name as
#' @family SpatialNetwork
#' @rdname SpatialNetwork-methods
setAs("SpatialNetwork", "SpatialPointsDataFrame",
	function(from) {
		x = V(from@g)$x
		y = V(from@g)$y
		SpatialPointsDataFrame(SpatialPoints(cbind(x,y), from@proj4string),
			data.frame(vertex.attributes(from@g)))
	}
)

#' As("SpatialNetwork", "igraph")
#'
#' @aliases coerce,SpatialNetwork,igraph-method
#' @name as
#' @family SpatialNetwork
#' @rdname SpatialNetwork-methods
setAs("SpatialNetwork", "igraph", function(from) from@g)

#' plot method for SpatialNetwork
#'
#' @name plot
#' @aliases plot,SpatialNetwork,missing-method
#' @docType methods
#' @param y ignored
#' @param arrow_size numeric; scaling factor for the arrow head
#' @rdname SpatialNetwork-methods
#' @examples
#' data(torontocentre)
#' tc = SpatialNetwork(torontocentre, direction = torontocentre$DIRECTION_)
#' library(sp)
#' plot(tc, arrow_size = .5)
#' points(tc)
setMethod("plot", signature(x = "SpatialNetwork", y = "missing"),
    function(x, y, ..., arrow_size = 0) {
		plot(as(x, "SpatialLines"), ...)
		if (arrow_size > 0)
			plotSpatialArrows(x, ..., arrow_size = arrow_size, direction = x$direction)
    }
)

## simplified from sp::plotSpatialLines
plotSpatialArrows <- function(SL, col = 1, lwd = 1, lty=1, 
	lend = 0, ljoin = 0, lmitre = 10, ..., arrow_size, direction) {

	n = length(SL@lines)
	if (length(col) != n) col <- rep(col, n)
	if (length(lwd) != n) lwd <- rep(lwd, n)
	if (length(lty) != n) lty <- rep(lty, n)

	for (i in seq(along = SL@lines)) { # grouping by link dir might speed this up.
		if (direction[i] != 0) {
			crds = SL@lines[[i]]@Lines[[1]]@coords
			if (direction[i] == -1) # up-link
				pt = crds[1:2,]
			if (direction[i] == 1)  # down-link
				pt = crds[nrow(crds):(nrow(crds)-1),]
     		arrows(pt[2,1], pt[2,2], pt[1,1], pt[1,2], length = 0.25 * arrow_size, angle = 25,
	             code = 2, col = col[i], lty = lty[i], lwd = lwd[i]) # ignore ...
		}
	}
}
