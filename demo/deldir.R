library(spnetwork)
library(deldir)
# The following function converts a set of points into a SpatialPolygons or into a SpatialLines object:

dd <- function(x, ..., to = "polygons") {
    stopifnot(is(x, "Spatial"))
    cc = coordinates(x)
    dd = deldir(list(x = cc[, 1], y = cc[, 2]), ...)
    if (to == "polygons") {
        tm = triMat(dd)
        fn = function(ix) {
            pts = tm[ix, ]
            pts = c(pts, pts[1])
            Polygons(list(Polygon(rbind(cc[pts, ]))), ix)
        }
        SpatialPolygons(lapply(1:nrow(tm), fn), proj4string = x@proj4string)
    } else if (to == "lines") {
        segs = dd$delsgs
        lst = vector("list", nrow(segs))
        for (i in 1:nrow(segs)) lst[i] = Lines(list(Line(cc[c(segs[i, 5], segs[i, 
            6]), ])), i)
        SpatialLines(lst, proj4string = x@proj4string)
    } else stop("argument to should be polygons or lines")
}
# We'll generate a set of 100 points in a unit square:

set.seed(5432)
x = runif(100)
x = x[order(x)]
y = runif(100)
library(sp)
pts = SpatialPoints(cbind(x, y))
# Next, we'll get the SpatialLines object from it:

sl = dd(pts, to = "lines")
## 
##      PLEASE NOTE:  The components "delsgs" and "summary" of the 
##      object returned by deldir() are now DATA FRAMES rather than 
##      matrices (as they were prior to release 0.0-18). 
##      See help("deldir").
##  
##      PLEASE NOTE: The process that deldir() uses for determining
##      duplicated points has changed from that used in version
##      0.0-9 of this package (and previously). See help("deldir").

# and plot it:

plot(sl)

# From this object, we can create a SpatialNetwork by

ln = SpatialNetwork(sl)
# when plotted, we can again use colour to denote the number of vertices connected to each edge:

plot(ln@g$x, ln@g$y, col = ln@g$n, pch = 16, cex = 2, asp = 1)
lines(sl, col = "grey")

# Now we compute the shortest path from the left-most point (1) to the right-most one (100):

library(igraph)
sp = as.vector(get.shortest.paths(ln@g, 1, 100)$vpath[[1]])
# and plot it

plot(ln@g$x, ln@g$y, col = ln@g$n, pch = 16, cex = 1.5, asp = 1)
lines(sl, col = "grey")
points(ln@g$x[sp], ln@g$y[sp], col = "red", cex = 2)
text(ln@g$x[c(1, 100)], ln@g$y[c(1, 100)], c("start", "end"), pos = 4)

# As the edge weights are computed by Line lengths, this is the geographically shortest path.
