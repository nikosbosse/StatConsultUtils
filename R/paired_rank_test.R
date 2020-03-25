#' Asymptomtic Paired Rank Test for two paired samples
#'
#' @description asymptotic rank test for two dependent samples that is based on
#' overall mid-ranks.
#'
#' Null hypothesis: The two samples are exchangeble, i.e. (X_k1, X_k2) is
#' equal in distribution to (X_k2, X_k1) where (X_k1, X_k2)' denotes independent
#' and identically distributed pairs of observations from k = 1,2,...,n
#' subjects.
#'
#' Applicability:
#' - Applicable for metric als well as ordinal data
#' - should only be used for n > 15
#'
#' Serves as a replacement for the paired t-test and the Wilcoxon signed rank
#' test that both can only be applied to metric data.
#'
#' Advantages of the paired rank test:
#' - invariant to monotonic transformations. Data e.g. can be logged or squared
#' without changing the results.
#'
#' Disadvantages:
#' - Less efficient than the t-test if it is applicable
#'
#' @param obs1 vector with observations of size k, where k is the number
#' of dependent samples.
#' @param obs2 vector of observations corresponding to obs1
#' @return A list
#' @importFrom stats t.test
#' @export
#' @examples
#'
#' ## Beispiel aus Munzel and Brunner, An Exact Paired Rank Test, 2002
#' obs1 <- c(6, 3, 5, 4, 5, 3, 4, 5, 5, 4, 6, 4, 4, 5, 6)
#' obs2 <- c(4, 1, 3, 4, 6, 3, 3, 4, 3, 3, 5, 5, 3, 4, 5)
#'
#' asymptotic_paired_rank_test(obs1, obs2)
#'

asymptotic_paired_rank_test <- function(obs1, obs2) {

  n <- length(obs1)

  # Gesamtmittelränge bestimmen und Rangdifferenzen
  rank <- rank(c(obs1, obs2))
  ranks_obs1 <- rank[1:n]
  ranks_obs2 <- rank[(n+1):(2 * n)]

  # t-Test auf den Rängen durchführen
  return(stats::t.test(ranks_obs1, ranks_obs2))

}








#' Recursive Paired Rank Test for two paired samples
#'
#' @description asymptotic rank test for two dependent samples that is based on
#' overall mid-ranks.
#'
#' Null hypothesis: The two samples are exchangeble, i.e. (X_k1, X_k2) is
#' equal in distribution to (X_k2, X_k1) where (X_k1, X_k2)' denotes independent
#' and identically distributed pairs of observations from k = 1,2,...,n
#' subjects.
#'
#' Applicability:
#' - Applicable for metric als well as ordinal data
#'
#' Serves as a replacement for the paired t-test and the Wilcoxon signed rank
#' test that both can only be applied to metric data.
#'
#' Advantages of the paired rank test:
#' - invariant to monotonic transformations. Data e.g. can be logged or squared
#' without changing the results.
#'
#' Disadvantages:
#' - Less efficient than the t-test if it is applicable
#'
#' @param obs1 vector with observations of size k, where k is the number
#' of dependent samples.
#' @param obs2 vector of observations corresponding to obs1
#' @return A list
#' @export
#' @examples
#'
#' ## Beispiel aus Munzel and Brunner, An Exact Paired Rank Test, 2002
#' obs1 <- c(6, 3, 5, 4, 5, 3, 4, 5, 5, 4, 6, 4, 4, 5, 6)
#' obs2 <- c(4, 1, 3, 4, 6, 3, 3, 4, 3, 3, 5, 5, 3, 4, 5)
#'
#' exact_paired_rank_recursive(obs1, obs2)
#'

## Rekursive Version des Tests
exact_paired_rank_recursive <- function(obs1, obs2) {
  n <- length(obs1)

  # Gesamtmittelränge bestimmen und Rangdifferenzen
  rank <- rank(c(obs1, obs2))
  ranks_obs1 <- rank[1:n]
  ranks_obs2 <- rank[(n+1):(2 * n)]
  D <- 2 * (ranks_obs2 - ranks_obs1)

  df <- data.frame(obs1, obs2, ranks_obs1, ranks_obs2, D)
  Tr <- sum(df$D)

  # Summe der positiven und negativen Rangdifferenzen
  Dplus <- sum(df$D[df$D > 0])
  Dminus <- sum(df$D[df$D < 0])

  # Funktion, um rekursiv p_n(t) zu bestimmen.
  p <- function(t, n, D) {
    if (n == 1 && t == 0) {
      return(0.5)
    } else if (n == 1 && t == abs(D[1])) {
      return(0.5)
    } else if (n == 1) {
      return(0)
    } else {
      res <- 0.5 *
        (p(n = n - 1,
           D = D,
           t = (t - abs(D[n]))) +
           p(n = n - 1,
             D = D,
             t = t))
      return(res)
    }
  }

  # Funktionswert p_n(t) bestimmen für alle t = 1:Dplus, mit n = Zahl der
  # verbundenen Beobachtungspaare, dann aufsumieren und mal zwei nehmen.
  pval <- 2 * sum(sapply(1:Dplus, FUN = p,
                         n = n, D = df$D))
  return(list(pval = pval,
              Dplus = Dplus,
              Tr = Tr))
}





#' Shift version Paired Rank Test for two paired samples
#'
#' @description Exact rank test for two dependent samples that is based on
#' overall mid-ranks. Applicable for metric als well as ordinal data
#'
#' Serves as a replacement for the paired t-test and the Wilcoxon signed rank
#' test that both can only be applied to metric data.
#'
#' Advantages of the paired rank test:
#' - invariant to monotonic transformations. Data e.g. can be logged or squared
#' without changing the results.
#'
#' @param obs1 vector with observations of size k, where k is the number
#' of dependent samples.
#' @param obs2 vector of observations corresponding to obs1
#' @return A list
#' @export
#' @importFrom dplyr lag
#' @examples
#'
#' ## Beispiel aus Munzel and Brunner, An Exact Paired Rank Test, 2002
#' obs1 <- c(6, 3, 5, 4, 5, 3, 4, 5, 5, 4, 6, 4, 4, 5, 6)
#' obs2 <- c(4, 1, 3, 4, 6, 3, 3, 4, 3, 3, 5, 5, 3, 4, 5)
#'
#' exact_paired_rank_shift(obs1, obs2)
#'


exact_paired_rank_shift <- function(obs1, obs2) {
  n <- length(obs1)

  # Gesamtmittelränge bestimmen und Rangdifferenzen
  rank <- rank(c(obs1, obs2))
  ranks_obs1 <- rank[1:n]
  ranks_obs2 <- rank[(n+1):(2 * n)]
  D <- 2 * (ranks_obs2 - ranks_obs1)

  df <- data.frame(obs1, obs2, ranks_obs1, ranks_obs2, D)
  Tr <- sum(df$D)

  # Summe der positiven und negativen Rangdifferenzen
  Dplus <- sum(df$D[df$D > 0])
  Dminus <- sum(df$D[df$D < 0])

  ## function to find greatest common divisor gcd
  ## größter gemeinsamer Teiler ggt
  gcd <- function(x) {
    k <- min(x[x != 0])
    while (k >= 1) {
      if (all(x/k == as.integer(x/k))) {
        return(k)
      } else {
        k <- k - 1
      }
    }
    return(1)
  }

  # Absolutwerte der Rangdifferenzen
  ad <- abs(df$D)
  ggt <- gcd(ad)

  # Teilne Vektor durch ggt um Rechenprozedur zu vereinfachen
  adshort <- ad / ggt

  # Größtmöglicher Wert, der erreicht werden kann
  lng <- sum(adshort)

  ## Shiftalgorithmus.
  # Starte mit Vektor c(1, 0, ..., 0)
  v <- c(1, rep(0, lng))

  # shifte Vektor um einen Wert = abs(adshort)[i] und addiere auf
  for (value in adshort) {
    v <- v + dplyr::lag(v, value, default = 0)
  }

  # relative Häufigkeiten bestimmen
  probs <- v / sum(v)
  vals <- ggt * 0:lng

  # Aufsummieren der relativen Häufigkeiten aller möglichen Werte
  # kleiner oder gleich D+
  pval <- 2 * sum(probs[vals <= Dplus])

  return(list(pval = pval,
              probs = probs,
              vals = vals,
              Dplus = Dplus,
              Tr = Tr))
}

