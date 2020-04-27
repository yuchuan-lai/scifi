\name{select}
\alias{select}

\title{
Select annual climate variable from the annual dataset
}
\description{
Processing downloaded annual climate records for a city for a particular temperature and precipitation variable.
}
\usage{
select(city.hist.data, var.name, year.col.name = "Year", thresholds = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{city.hist.data}{
A dataframe object that contains the different annual temperature and precipitaiton variable and years for a city. It can be a datafram object that is downloaded using \code{\link{download}} function.
}
  \item{var.name}{
A character object that specifies the column name of the selected climate variable in the city.hist.data.
}
  \item{year.col.name}{
A character object that specifies the column name used for the information of years in the city.hist.data.
}
  \item{thresholds}{
A number specifies the number of missing daily data acceptable for the selection of annual climate variable.
}
}

\value{
A dataframe object that provides the information about the years and one selected climate variable for the city.
}
\author{
Yuchuan Lai
}

\seealso{
See Also as \code{\link{citylist}} and \code{\link{download}}.
}
\examples{
# Download the annual data for Pittsburgh
pit.annual <- download("Pittsburgh", "annual")

# Process the annual data of Pittsburgh to obtain annual average temperature as a data.frame object
pit.temp <- select(pit.annual, "Avg.Temp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
