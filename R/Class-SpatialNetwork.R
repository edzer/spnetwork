#' @import methods sp igraph
#' @importFrom stats weights
#' @importFrom graphics arrows points

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
#'    \item{\code{weightfield}}{character; describing the weight field used}
#'  }
#'
#' @usage SpatialNetwork(sl, g, weights, weightfield, directions, zero)
#' @param sl object of one of (a sublasses of) \link{SpatialLines}, with links
#' @param g object of class \link{igraph}; if missing, the graph is sorted out from the links
#' @param weights weight for links (defaults to length of linestring)
#' @param weightfield character; name of the attribute field of \code{sl} that will be used as weights
#' @param directions numeric; if omitted, undirected graph, else integer vector indicating direction of a link: -1 for up-link, 0 for two-way, or 1 for down-link
#' @param zero numeric; zero threshold passed on to \link[sp]{zerodist}
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
#' plot(sln)
#' points(sln)
#' library(igraph) # E
#' text(V(sln@g)$x, V(sln@g)$y, E(sln@g), pos = 4)
#' plot(sln@g)
setClass("SpatialNetwork",
	contains = "SpatialLinesDataFrame", 
	slots = c(g = "igraph", weightfield = "character"),
    validity = function(object) {
		# print("Entering validation: SpatialNetwork")
		if (any(sapply(object@lines, function(x) length(x@Lines)) != 1))
			stop("all Lines objects need to have a single Line") 
    	if (!is_directed(object@g)) {
			if (length(object@lines) != length(E(object@g)))
				stop("edges do not match number line segments") 
			if (!identical(E(object@g)$link_index, 1:length(object@lines))) {
				print(E(object@g)$link_index)
				stop("link_index should equal 1:length(lines)") 
			}
		} else { # directed:
			if (length(object@lines) > length(E(object@g)))
				stop("# edges should not not be smaller than the number of line segments") 
		}
		return(TRUE)
	}
)

#' @export
SpatialNetwork = function(sl, g, weights, weightfield, directions, zero = 0.0) {
    stopifnot(is(sl, "SpatialLines"))
    if (!is(sl, "SpatialLinesDataFrame")) 
        sl = new("SpatialLinesDataFrame", sl, data = data.frame(id = 1:length(sl)))
	if (missing(g)) { # sort out from sl
    	startEndPoints = function(x) {
        	firstLast = function(L) {
            	cc = coordinates(L)[[1]]
            	rbind(cc[1, ], cc[nrow(cc), ])
        	}
        	do.call(rbind, lapply(x@lines, firstLast))
    	}
    	s = startEndPoints(sl)
    	zd = zerodist(SpatialPoints(s, sl@proj4string), zero)
    	pts = 1:nrow(s) # 1-2 3-4 5-6 etc
	
		# replace higher with lower, identical points, e.g. 1-2 3-1 2-6
    	# the following can't be done vector-wise, there is a progressive effect:
    	if (nrow(zd) > 0) {
        	for (i in 1:nrow(zd)) 
				pts[zd[i, 2]] = pts[zd[i, 1]]
    	}
    	stopifnot(identical(s, s[pts, ]))
	
    	# map to 1:length(unique(pts)), e.g. 1-2 3-1 2-4
    	pts0 = match(pts, unique(pts))
		link_index = 1:length(sl)
		if (! missing(directions)) { # -1: upstream, 0: two-way, 1: down-stream
			stopifnot(length(directions) == length(sl))
			stopifnot(all(directions %in% c(-1,0,1)))
			sl$directions = directions
			downlink = directions != -1 # 1..length(directions)
			e = matrix(pts0, ncol = 2, byrow = TRUE) # each row is an edge/link
			# handle -1: reversed edges for one-way, up-link
			rev = which(directions == -1)
			to = e[rev,2]
			e[rev,2] = e[rev,1]
			e[rev,1] = to
			# 0: add reverse direction edges
			two_way = which(directions == 0)
			if (length(two_way) > 0) {
				e = rbind(e, e[two_way, c(2,1)]) # add up-links
				downlink = c(downlink, rep(FALSE, length(two_way)))
				link_index = c(link_index, two_way)
			}
			# 1: all is fine, do nothing (one-way, down-link)
    		g = graph(t(e), directed = TRUE)  # edges
        	E(g)$downlink = downlink
		} else
    		g = graph(pts0, directed = FALSE)  # edges
		E(g)$link_index = link_index
    	nodes = s[unique(pts), ]
		# this needs work:
    	V(g)$x = nodes[, 1]  # x-coordinate vertex
    	V(g)$y = nodes[, 2]  # y-coordinate vertex
    	V(g)$n = as.vector(table(pts0))  # nr of edges
    	# line lengths:
    	if (missing(weights) && missing(weightfield)) {
    		sl$length = SpatialLinesLengths(sl) # takes care of projected/geographical
    	    weightfield = 'length'
		}
    	# create list with vertices, starting/stopping for each edge?  add for
    	# each SpatialLines, the start and stop vertex
    	pts2 = matrix(pts0, ncol = 2, byrow = TRUE)
    	sl$start = pts2[, 1]
    	sl$end = pts2[, 2]
	}
	if (! missing(weights)) {
		stopifnot(length(weights) == length(sl))
	    if (missing(weightfield))
	        weightfield = 'weight'
        sl@data[, weightfield] <- weights
	} else
        weights = sl@data[, weightfield]
    E(g)$weight = weights[E(g)$link_index]
    new("SpatialNetwork", sl, g = g, weightfield = weightfield)
}
