\name{OLYMPIC_DECATHLON_2016}
\alias{OLYMPIC_DECATHLON_2016}
\docType{data}
\title{OLYMPIC_DECATHLON_2016
}
\description{Results of the athletes competing in the decathlon at the Games of the XXXI Olympiad (Rio de Janeiro, Brazil, year 2016). The dataset contains the points scored in each event by the 23 decathletes who who completed the dacathlons, together with their names and nationalities. It is freely available at www.iaaf.org, the official website of the IAAF (International Association of Athletics Federations).
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("OLYMPIC_DECATHLON_2016")}
\format{
  A data frame with 23 observations on the following 13 variables.
  \describe{
    \item{\code{OS}}{a numeric vector. Athletes' ranking.}
    \item{\code{ATHLETE}}{a factor with levels \code{Adam SebastianHELCELET} \code{Akihiko NAKAMURA} \code{Arthur ABELE} \code{Ashton EATON} \code{Bastien AUZEIL} \code{Cedric DUBLER} \code{Damian WARNER} \code{Dominik DISTELBERGER} \code{Jeremy TAIWO} \code{Kai KAZMIREK} \code{Karl Robert SALURI} \code{Keisuke USHIRO} \code{Kevin MAYER} \code{Kurt FELIX} \code{Larbi BOURRADA} \code{Leonel SUAREZ} \code{Lindon VICTOR} \code{Luiz Alberto DE ARAUJO} \code{Pau TONNESEN} \code{Pawel WIESIOLEK} \code{Thomas VAN DER PLAETSEN} \code{Yordani GARCIA} \code{Zach ZIEMEK}}
    \item{\code{COUNTRY}}{a factor with levels \code{ALG} \code{AUS} \code{AUT} \code{BEL} \code{BRA} \code{CAN} \code{CUB} \code{CZE} \code{ESP} \code{EST} \code{FRA} \code{GER} \code{GRN} \code{JPN} \code{POL} \code{USA}}
    \item{\samp{100.METRES}}{a numeric vector. Points scored in the one hundred metres.}
    \item{\code{LONG.JUMP}}{a numeric vector. Points scored in the long jump.}
    \item{\code{SHOT.PUT}}{a numeric vector. Points scored in the shot put.}
    \item{\code{HIGH.JUMP}}{a numeric vector. Points scored in the high jump.}
    \item{\samp{400.METRES}}{a numeric vector. Points scored in the four hundred metres.}
    \item{\samp{110.METRES.HURDLES}}{a numeric vector. Points scored in the one hundred and ten metres hurdles.}
    \item{\code{DISCUS.THROW}}{a numeric vector. Points scored in the discus throw.}
    \item{\code{POLE.VAULT}}{a numeric vector. Points scored in the pole vault.}
    \item{\code{JAVELIN.THROW}}{a numeric vector. Points scored in the javelin throw.}
    \item{\code{X1500.METRES}}{a numeric vector. Points scored in the one thousand and five hundred metres.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{www.iaaf.org
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(OLYMPIC_DECATHLON_2016)
## maybe str(OLYMPIC_DECATHLON_2016) ; plot(OLYMPIC_DECATHLON_2016) ...
}
\keyword{datasets}
