#' Control Parameters for QUINT Algorithm
#'
#' Various parameters that control aspects of the \dQuote{quint} algorithm.
#' Appendix A of Dusseldorp & Van Mechelen (2013) gives a detailed overview
#' of the choices that can be made.
#'
#' @param crit the type of difference in treatment outcome used in the
#'   partitioning criterion: "es" (Treatment effect sizes) or
#'   "dm" (Difference in treatment means).
#' @param maxl maximum number of leaves (\eqn{L}) of the tree. Default value is 10.
#' @param a1 the minimal sample size of Treatment A (\eqn{T}=1) in a leaf.
#' If NULL, a1 is set to 1/10 of the sample size of the Treatment A group (assignment is done in the function quint). The minimum value is 2.
#' @param a2 the minimal sample size of Treatment B (\eqn{T}=2) in a leaf.
#' If NULL, a2 is set to 1/10 of the sample size of the Treatment B group (assignment is done in the function quint). The minimum value is 2.
#' @param w a vector with w1 and w2 representing the weights of, respectively,
#'   the Difference in treatment outcome component and the Cardinality component
#'   of the partitioning criterion. If crit = "dm", the default value of \eqn{w1}
#'   is 1/ \eqn{log}(1+IQR(Y)). If crit = "es", the default value of \eqn{w1} is
#'   1/\eqn{log}(1+3). The default of \eqn{w2} is 1/\eqn{log}(0.50\eqn{N}).
#' @param Bootstrap whether the bias-corrected bootstrap procedure should
#'   be performed. The default is TRUE.
#' @param B the number of bootstrap samples to be drawn. The default is 25. We recommend a number of bootstraps of at least 25.
#' @param dmin the minimum absolute standardized mean difference in
#'  treatment outcome in one of the leaves assigned to P1 and one of the leaves
#'  assigned to P2 of the pruned tree. This value is used to check whether a qualitative interaction
#'  is present in the data (the qualitative interaction condition); dmin controls
#'  the balance between Type I error and Type II error. The default value of dmin is 0.30.
#'
#' @return A list containing the options.
#' @references Dusseldorp, E., Doove, L., & Van Mechelen, I. (2016). Quint:
#'   An R package for the identification of subgroups of clients who differ in
#'   which treatment alternative is best for them. \emph{Behavior Research Methods,
#'   48}(2), 650-663. DOI 10.3758/s13428-015-0594-z
#'
#'   Dusseldorp E. and Van Mechelen I. (2014). Qualitative interaction
#'   trees: a tool to identify qualitative treatment-subgroup interactions.
#'   \emph{Statistics in Medicine, 33}(2), 219-237. DOI: 10.1002/sim.5933.
#'
#' @seealso \code{\link{quint}}
#' @examples data(bcrp)
#' formula1<- I(cesdt1-cesdt3)~cond | nationality+marital+wcht1+age+
#'   trext+comorbid+disopt1+uncomt1+negsoct1
#' #Specify the Difference in treatment outcome as Difference in means
#' #and skip the bias-corrected bootstrap procedure
#' #and change the maximum number of leaves
#' control3<-quint.control(crit="dm",Bootstrap=FALSE,maxl=3)
#' quint3<-quint(formula1, data= subset(bcrp,cond<3),control=control3)
#' summary(quint3)
#'
#' #Set number of bootstrap samples at 30
#' control4<-quint.control(B=30)
#'
#' #Set minimal sample size in each treatment group at 5
#' control5<-quint.control(a1=5,a2=5)
#'
#' @export
quint.control<-function(crit="es", maxl = 10, a1=NULL, a2=NULL,
                         w=NULL, Bootstrap=TRUE, B=25, dmin=0.30) {
  if(crit != 'es' & crit!= 'dm'){
    stop("Criterion ", crit, " specified is not effect size (es) or difference in means (dm).")}
  #crit="es" (effect size criterion) or "dm" (difference in means criterion)

  if(!is.null(a1)){
    parvec<-round(c(a1,a2))
    if(a1<2){warning("a1 should be grater or equal than 2.")}
      if(!is.null(a2)){
        if(a2<2){warning("a2 should be grater or equal than 2.")}
      }
  }else{
    if(!is.null(a2)){
      parvec<-round(c(a1,a2))
      if(a2<2){warning("a2 should be grater or equal than 2.")}
    }else{
      parvec<-NULL
    }
  }



  list(crit = crit, maxl = maxl, parvec = parvec, w = w, Boot = Bootstrap, B = B, dmin = dmin)
}
