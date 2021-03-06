% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_chromatogram.R
\name{addFlankToLeft}
\alias{addFlankToLeft}
\title{Add signal to the left of XIC}
\usage{
addFlankToLeft(flankSeq, XIC, newXIC)
}
\arguments{
\item{flankSeq}{(logical) must be TRUE at the front of the vector.}

\item{XIC}{(data-frame) first column is time, second column is intensity.}

\item{newXIC}{(data-frame) first column is time, second column is intensity.}
}
\value{
(dataframe) has two columns:
\item{time}{(numeric)}
\item{intensity}{(numeric)}
}
\description{
This function copies flanking chromatogram from XIC and paste it to the left of newXIC.
}
\examples{
time <- seq(from = 3003.4, to = 3048, by = 3.4)
y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
       4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
chrom <- data.frame(time, y)
chrom2 <- data.frame(time = c(3013.4, 3016, 3020), intensity = c(1.2, 3.4, 5.6))
flankSeq <- as.logical(c(1,1,0,0,0,0,0,0,0,0,0,0,1,1))
\dontrun{
addFlankToLeft(flankSeq, chrom, chrom2)
}
}
\seealso{
\code{\link{childXIC}, \link{addFlankToRight}}
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2020) + GPL-3
Date: 2020-07-16
}
\keyword{internal}
