\name{tumaGUI}
\alias{tumaGUI}
\title{Start the Tuma GUI}
\usage{
  tumaGUI(data = borenstein_tuma)
}
\arguments{
  \item{data}{an optional data frame. If a data frame is
  provided, it has to have at least the following three
  elements: T, se, studlab. If the name of a newly created
  data frame ends with "_tuma" (e.g., "borenstein_tuma",
  "df_tuma", "mydata_tuma" etc.), it will show up in the
  drop down menu.}
}
\description{
  The Tuma GUI can be started by typing \code{tumaGUI()}.
  In its current version, the package comes with two data
  sets (data.frame objects): \itemize{ \item
  borenstein_tuma (Source: Borenstein, M., Hedges, L. V.,
  Higgins, J. P. T., & Rothstein, H. R. (2009).
  Introduction to meta-analysis. John Wiley and Sons.
  (Table 14.7, p. 98) \item raudenbush_tuma (Source:
  Raudenbush, S. W. (1984). Magnitude of Teacher Expectancy
  Effects on Pupil IQ as a Function of the Credibility of
  Expectancy Induction: A Synthesis of Findings From 18
  Experiments. Journal of Educational Psychology, 76(1),
  85--97. }
}
\examples{
tumaGUI()
tumaGUI(raudenbush_tuma)
}
\author{
  Bernd Weiss \email{bernd.weiss@uni-koeln.de}
}

