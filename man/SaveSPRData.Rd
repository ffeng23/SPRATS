% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPRTwoState.R
\name{SaveSPRData}
\alias{SaveSPRData}
\title{S3 function to save SensorgramData data}
\usage{
SaveSPRData(DATA, file, sep = "\\t")
}
\arguments{
\item{DATA}{SensorgramData to be saved}

\item{file}{the directory path to the file}

\item{sep}{character to delimit the text fields}
}
\description{
save SensorgramData object data to the disk.
}
\details{
It save the object data to the disk in a text format,
	following the SensiQ format. See the detail \code{\link{ReadSensorgramData}}
	Note: it will issue a warning to complain about appending the table.
	There is nothing wrong with it. Simple ignore this warning.
}
\seealso{
\code{\link{ReadSensorgramData}} \code{\link{SensorgramData-class}} \code{\link{plot}}
}

