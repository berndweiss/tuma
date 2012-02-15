\name{initTclVarLists}
\alias{initTclVarLists}
\title{Initialize lists of tcl variables containing effect sizes, standard errors
and study labels}
\usage{
  initTclVarLists(T, se, studlab)
}
\arguments{
  \item{T}{numeric vector of effect sizes}

  \item{se}{numeric vector of standard errors}

  \item{studlab}{character vector of study labels}
}
\description{
  This function initializes several lists of tcl variables.
  Note that the assignment has to take place in the global
  environment (i.e. use "<<-"). In particular, the
  following lists of tcl variables are generated: \itemize{
  \item lstTclVarT \item lstTclVarSe \item lstTclVarStudlab
  \item lstTclEnT \item lstTclEnSe \item lstTclEnStudlab
  \item lstTclSlT \item lstTclSlSe }
}
\keyword{internal}

