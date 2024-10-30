#' @importFrom formula.tools env
#' @importFrom stats model.frame
#' @importFrom stats terms
#' @importFrom stats model.response
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
formula.to.data <- function(formula,
                            data, 
                            simulations   = FALSE,
                            fixed.effects = FALSE) {
  
  ## Extract data from the formula
  if (missing(data)) {
    data      <- env(formula)
  }
  formula     <- as.formula(formula)
  yname       <- NULL
  
  if (simulations) {
    if(length(formula) != 2) stop("formula is not correct. To simulate data the expected format is ~ X1 + X2 + ....")
  } else {
    if(length(formula) != 3) stop("formula is not correct. For the estimation, the expected format is y ~ X1 + X2 + ....")
    yname     <- all.vars(formula)[1]
  }
  
  ## call model.frame()
  mf          <- model.frame(formula, data = data)
  ## extract response, terms, model matrices
  y           <- model.response(mf, "numeric")
  X           <- model.matrix(terms(formula, data = data, rhs = 1), mf)
  xname       <- colnames(X)
  intercept   <- "(Intercept)" %in% xname
  if(fixed.effects & intercept){
    X         <- X[, xname != "(Intercept)", drop = FALSE]
    xname    <- xname[xname != "(Intercept)"]
    intercept <- FALSE
  }
  
  list("formula"   = formula, 
       "X"         = X, 
       "y"         = y,
       "intercept" = intercept,
       "yname"     = yname,
       "xname"     = xname)
}
