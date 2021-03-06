#' R package for Qualitative Treatment-Subgroup Interactions
#'
#' When two treatment alternatives (say A and B) are available  for some problem,
#' one may be interested in qualitative treatment-subgroup interactions. Such
#' interactions imply the existence of subgroups of persons (patients) which are
#' such that in one subgroup Treatment A outperforms Treatment B, whereas the reverse
#' holds in another subgroup. Obviously, this type of interactions is crucial for
#' optimal treatment assignment of future patients. Given baseline characteristics and
#' outcome data from a two-arm Randomized Controlled Trial (RCT), QUalitative INteraction
#' Trees (QUINT) is a tool to identify subgroups that are involved in meaningful
#' qualitative treatment-subgroup interactions. The result of QUINT is a tree that
#' partitions the total group of participants (patients) on the basis of their baseline
#' characteristics into three subgroups (i.e., partition classes): Subgroup 1: Those for
#' whom Treatment A is better than Treatment B (P1), Subgroup 2: Those for whom
#' Treatment B is better than Treatment A (P2), and Subgroup 3: Those for whom it does
#' not make any difference (P3).
#'
#' @details \tabular{ll}{
#'   Package: \tab quint\cr
#'   Type: \tab Package\cr
#'   Version: \tab 2.2.0\cr
#'   Date: \tab 2020-02-03\cr
#'   License: \tab GPL\cr
#' }
#'
#'   This method is suitable for a continuous outcome variable. From version 1.2 onwards
#'   the baseline variables for growing a tree may have numerical
#'   or integer values (such as continuous, ordinal or dichotomous variables) or may be nominal
#'   (categorical variables with factors). Previously only numerical or dichotomous variables
#'   were supported. Another new feature of this version is that
#'   the output of a \code{quint} object can now also display results for either the raw difference
#'   in means or the effect size with corresponding standard error. This depends on the criterion
#'   specified. Furthermore a predict function \code{predict.quint} is newly included in this
#'   package. The final new feature is a validate function \code{quint.validate} for estimating
#'   the bias (i.e., optimism) of a grown QUINT tree.
#'
#'   From version 2.0 onwards the qualitative treatment-subgroup interaction is checked during the prune
#'   of the tree and not at the begining of QUINT. Furthermore, it is possible to obtain outcomes
#'   from the summary and predict functions when the tree only contains the root node.
#'
#'   The core function of the package is \code{\link{quint}}.
#'
#' @author Maintainer: Elise Dusseldorp <elise.dusseldorp@fsw.leidenuniv.nl>
#' @references Dusseldorp, E., Doove, L., & Van Mechelen, I. (2016). Quint:
#'   An R package for the identification of subgroups of clients who differ in
#'   which treatment alternative is best for them. \emph{Behavior Research Methods,
#'   48}(2), 650-663. DOI 10.3758/s13428-015-0594-z
#'
#'   Dusseldorp E. and Van Mechelen I. (2014). Qualitative interaction
#'   trees: a tool to identify qualitative treatment-subgroup interactions.
#'   \emph{Statistics in Medicine, 33}(2), 219-237. DOI: 10.1002/sim.5933.
#'
#'   Scheier M.F., Helgeson V.S., Schulz R., et al.(2007). Moderators of interventions
#'   designed to enhance physical and psychological functioning among younger women
#'   with early-stage breast cancer. \emph{Journal of Clinical Oncology, 25}, 5710-5714.
#'   DOI: 10.1200/JCO.2007.11.7093.
#' @keywords package
#' @seealso \code{\link{quint}},\code{\link{summary.quint}},\code{\link{quint.control}},
#'   \code{\link{prune.quint}},\code{\link{predict.quint}},\code{\link{quint.validate}},
#'   \code{\link{quint.bootstrapCI}}
#'
#' @docType package
#' @name quint-package
NULL
