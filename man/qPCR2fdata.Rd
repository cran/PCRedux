% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qPCR2fdata.R
\name{qPCR2fdata}
\alias{qPCR2fdata}
\title{A helper function to convert amplification curve data to the fdata format.}
\usage{
qPCR2fdata(data, preprocess = FALSE)
}
\arguments{
\item{data}{is a data set containing the amplification cycles (1. column)
and the fluorescence (subsequent columns).}

\item{preprocess}{is a logical parameter (default FALSE). If TRUE, the \code{CPP}
function from the chipPCR package (Roediger et al. 2015) is used to pre-process
the data (e.g., imputation of missing values).
and the fluorescence (subsequent columns).}
}
\value{
gives an \code{fdata} object (S3 class, type of \code{list}) as output 
for a converted amplification curve.
}
\description{
\code{qPCR2fdata} is a helper function to convert qPCR data to the functional
\code{fdata} class as proposed by  Febrero-Bande & de la Fuente (2012). This
function prepares the data for further analysis with the \code{fda.usc} package,
which includes utilities for functional data analysis (e.g., Hausdorff
distance).
}
\examples{
library(qpcR)
library(fda.usc)

default.par <- par(no.readonly = TRUE)
# Calculate slope and intercept on noise (negative) amplification curve data
# for the last eight cycles.
# Convert the qPCR data set to the fdata format
res_fdata <- qPCR2fdata(testdat)

# Extract column names and create rainbow color to label the data
res_fdata_colnames <- colnames(testdat[-1])
data_colors <- rainbow(length(res_fdata_colnames), alpha=0.5)

# Plot the converted qPCR data
par(mfrow=c(1,2))
plot(res_fdata, xlab="cycles", ylab="RFU", main="testdat", type="l",
                   lty=1, lwd=2, col=data_colors)
legend("topleft", as.character(res_fdata_colnames), pch=19,
         col=data_colors, bty="n", ncol=2)

# Calculate the Hausdorff distance (fda.usc) package and plot the distances
# as clustered data.

res_fdata_hclust <- metric.hausdorff(res_fdata)
plot(hclust(as.dist(res_fdata_hclust)), main="Clusters of the amplification\n
   curves as calculated by the Hausdorff distance")
par(default.par)
}
\references{
M. Febrero-Bande, M.O. de la Fuente, others, \emph{Statistical
computing in functional data analysis: The R package fda.usc}, Journal of
Statistical Software. 51 (2012) 1--28. http://www.jstatsoft.org/v51/i04/

S. Roediger, M. Burdukiewicz, P. Schierack, \emph{chipPCR: an R package to
pre-process raw data of amplification curves}, Bioinformatics. 31 (2015)
2900--2902. doi:10.1093/bioinformatics/btv205.
}
\author{
Stefan Roediger, Michal Burdukiewcz
}
\keyword{fdata}
