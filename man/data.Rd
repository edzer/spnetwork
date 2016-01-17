\name{ottawacentre}
\alias{ottawacentre}
\alias{torontocentre}
\docType{data}
\title{Ontario city centre networks}
\description{Ontario city centre network sample data for Ottawa and Toronto, in long/lat WGS84}
\format{
  Objects of class \link[sp]{SpatialLinesDataFrame} that contain as attribution amonst others
  \describe{
    \item{DIRECTION_}{numeric; indicating link direction: 0 two-way, 1 one-way down-link, -1 one-way up-link}
  }
}
\usage{
data(ottawacentre)
data(torontocentre)
}
\note{ Distributed under the Ontario Open Data Licence
(https://www.ontario.ca/page/open-government-licence-ontario) }

\author{ Network selection done by Richard Ellison; imported from shapefiles using \link[rgdal]{readOGR} }

\keyword{datasets}

\examples{
library(sp)
data(ottawacentre)
oc = SpatialNetwork(ottawacentre, directions = ottawacentre$DIRECTION_)
plot(oc)
points(oc)
data(torontocentre)
tc = SpatialNetwork(torontocentre, directions = torontocentre$DIRECTION_)
plot(tc)
points(tc)
}
