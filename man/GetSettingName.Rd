% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation.R
\name{GetSettingName}
\alias{GetSettingName}
\title{Helper function for simulations}
\usage{
GetSettingName(settings, digits = 3, display = c("short", "pretty", "tiny"))
}
\arguments{
\item{settings}{A named list. The names corresponding to the setting parameter names, and the values are the parameter values}

\item{digits}{How many digits to use to format the numerical values}

\item{display}{Which version of the setting name to produce}
}
\value{
A character vector of setting names
}
\description{
Get the name of the settings from a named list good for saving into a data frame or use as file name
}
