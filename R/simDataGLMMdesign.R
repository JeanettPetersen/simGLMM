#' Simulated data to use in GLMM
#' 
#' Here, a design matrix is specified for a univariate model, whereas in the function simDataGLMM, a fixed mean is used, and a multivariate model can be specified.
#' 
#' For the Gaussian family, the variance parameter is 1. For the Gamma distribution, the rate parameter is 1.
#'
#' @param fixedMat The model matrix (numeric matrix or data.frame of factors and numeric variables) of the fixed effects to enter the linear predictor. Make sure that the contrasts fit together with the specification of the coefficients in fixedCoef.
#' @param fixedCoef The coefficients to use with the fixed effects.
#' @param randomMat The model matrix (data.frame of factors) of the random effects. 
#' @param covM Covariance matrix of the random effects. Currently, only the diagonal is used.
#' @param disFam The distribution family. The canonical link function is chosen.
#' @param detailed.output Output all the steps. Default: F. (This option is only for testing and can be removed.)
#' @param seed integer. Seed used for the simulation. Default: NULL, so the seed provided by the system is used. 
#' 
#' @examples
#' 
#' # one numeric explanatory variable and 1 factor with 2 levels
#' fixedMat = as.matrix(runif(10))
#' randomMat = data.frame(group = as.factor(rep(c("A", "B"), each=5)))
#' simDataGLMMdesign(fixedMat, fixedCoef=0.5, randomMat, covM = 1, disFam = poisson(), seed = 1234)
#' 
#' fixedMat = as.matrix(runif(20), ncol=2)
#' randomMat = as.factor(rep(c("A", "B"), each=5))
#' simDataGLMMdesign(fixedMat, fixedCoef=0.5, randomMat, covM = 1, disFam = poisson())
#' 
#' fixedMat = data.frame(quant=runif(10), fac=as.factor(rep(c("D1", "D2"), each=5)))
#' randomMat = data.frame(group = as.factor(rep(c("A", "B"), each=5)), group2 = as.factor(rep(c("a", "b", "c", "d", "e"), 2)))
#' simDataGLMMdesign(fixedMat, fixedCoef=c(1, 5, 10), randomMat, covM = diag(c(1,2)), disFam = inverse.gaussian(), detailed.output=T)
#' 
#' @return 
#' A data frame consisting of the fixed effects and random effects as given in the input and the simulated response variable 'res'.


simDataGLMMdesign <- function(fixedMat, fixedCoef, randomMat, covM, disFam, detailed.output=F, seed=NULL) {
  
  require(mvtnorm)
  require(statmod)
  
  # sample size
  N = nrow(fixedMat)
  
  # Make design matrices, if needed.

  # fixedMat: 
  


  if(is.data.frame(fixedMat)){
    
    fixedForm = formula( paste("~ 0 +", paste(names(fixedMat), collapse='+')))
    fixedMM = model.matrix(fixedForm, fixedMat)
    
  } else if(is.matrix(fixedMat)) {
    
     fixedMM = fixedMat

  } else {
    
    stop("Please provide a dataframe or a numeric matrix as fixedMat.")
    
  }

  
  # randomMat:
  
  randomMM = matrix(NA, ncol=0, nrow=N)

  if(is.data.frame(randomMat)){
    
    for(i in 1:ncol(randomMat)){
      
      randomForm = formula( paste("~ 0 +", names(randomMat[i])))
      randomMM = cbind(randomMM, model.matrix(randomForm, randomMat))     
      
    }
    
  } else {
    stop("Please provide a dataframe as randomMat.")
  }
  
  
  # simulate random effects
  
  if(!is.null(seed)) set.seed(seed)
  
  randomVec = numeric(0)
  for(i in 1:ncol(randomMat)){
    randomVec = c( randomVec, rnorm(nlevels(randomMat[,i]), mean=0, sd=sqrt(diag(covM)[i])) )
  }
  
  
  # compute linear Part
  
  randomPart = randomMM %*% randomVec
  fixedPart = fixedMM %*% fixedCoef
  linPred = randomPart + fixedPart
  
  # simulate response 
  
  if (disFam$family == "poisson") {
    sim <- rpois(N, exp(linPred))
  } else if (disFam$family == "gaussian") {
    sigma = 1
    sim <- rnorm(N, linPred, sigma)
  } else if (disFam$family == "binomial") { 
    sim <- rbinom(N, 1, exp(linPred)/(1 + exp(linPred)))
  } else if (disFam$family == "Gamma") {
    rate=1
    sim <- rgamma(N, shape = (linPred)^(-1)/rate,
                  rate = rate)
  } else if (disFam$family == "inverse.gaussian") {
    sim <- rinvgauss(N, linPred^(-1/2))
  } else {
    stop("Family not known")
  }
  
  
  # output
  
  if(detailed.output){
    return(list(fixedMM = fixedMM, randomMM = randomMM, randomVec = randomVec, linPred = linPred, res = sim))
  } else {  
    return(data.frame(fixedMat, randomMat, res=sim))
  }
  
}
