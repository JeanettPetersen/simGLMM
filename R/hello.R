# Some useful keyboard shortcuts for package authoring: Build
# and Reload Package: 'Cmd + Shift + B' Check Package: 'Cmd +
# Shift + E' Test Package: 'Cmd + Shift + T'

# Function
# ----------------------------------------------------------------

#' Simulated data to use in GLMM
#'
#' @param clus The number of clusters wanted.
#' @param rep The number of repetitions in each cluster.
#' @param var The number of variables. Only needed in case of multivariate data, default is 1.
#' @param fixedMean A fixed value to enter the linear predictor. In case of multivariate data, the dimension needs to equal the number of variables.
#' @param covM The variance of the random component representing variation between clusters. In case of multidimensional data, the dimension needs to equal the number of variables.
#' @param disFam The distribution family. The canonical link function is chosen.
#' @examples
#' # Simulates 3 blocks with 5 repetitions in each. The linear predictor is 10+random component with variance 5.
#' simDataGLMM(clus=3,rep=5,fixedMean=10,covM=5,disFam = poisson())
#' # Simulates 3 blocks with 5 repetitions in each for 3 variables. The linear predictor is 10+random component with defined covariance matrix.
#' A <- matrix(runif(3^2)*2-1, ncol=3)
#' sigma <- t(A) %*% A
#' simDataGLMM(clus=3,rep=5,fixedMean=rep(10,3),covM=sigma,disFam = poisson())


simDataGLMM <- function(clus, rep, var = 1, fixedMean, covM, disFam) {
    require(mvtnorm)
    require(statmod)
    # Simulate one random efffect per cluster per variable
    if (var == 1) {
        RandomEffects <- matrix(rnorm(clus, 0, covM), ncol = 1)
        res <- matrix(nrow = clus * rep, ncol = var)
    } else {
        RandomEffects <- rmvnorm(clus, rep(0, var), covM)
        res <- matrix(nrow = clus * rep, ncol = var)
    }


    if (disFam$family == "poisson") {
        # Simulate data in each cluster for each variable
        for (i in 1:clus) {
            sim <- rpois(rep * var, exp(fixedMean + RandomEffects[i,
                ]))
            for (j in 1:var) {
                res[(1:rep) + rep * (i - 1), j] <- sim[seq(1, (var *
                  rep), var) + (j - 1)]
            }
        }
    } else if (disFam$family == "gaussian") {
        # Simulate data in each cluster for each variable
        for (i in 1:clus) {
            sim <- rnorm(rep * var, fixedMean + RandomEffects[i,
                ], runif(1))
            for (j in 1:var) {
                res[(1:rep) + rep * (i - 1), j] <- sim[seq(1, (var *
                  rep), var) + (j - 1)]
            }
        }
    } else if (disFam$family == "binomial") {
        for (i in 1:clus) {
            linPred <- fixedMean + RandomEffects[i, ]
            # Simulated binary numbers
            sim <- rbinom(rep * var, 1, exp(linPred)/(1 + exp(linPred)))
            for (j in 1:var) {
                res[(1:rep) + rep * (i - 1), j] <- sim[seq(1, (var *
                  rep), var) + (j - 1)]
            }
        }
    } else if (disFam$family == "Gamma") {
        for (i in 1:clus) {
            linPred <- fixedMean + RandomEffects[i, ]
            rate <- 1
            # Simulated binary numbers
            sim <- rgamma(rep * var, shape = (linPred)^(-1)/rate,
                rate = rate)
            for (j in 1:var) {
                res[(1:rep) + rep * (i - 1), j] <- sim[seq(1, (var *
                  rep), var) + (j - 1)]
            }
        }
    } else if (disFam$family == "inverse.gaussian") {
        for (i in 1:clus) {
            linPred <- fixedMean + RandomEffects[i, ]
            # Simulated binary numbers
            sim <- rinvgauss(rep * var, linPred^(-1/2))
            for (j in 1:var) {
                res[(1:rep) + rep * (i - 1), j] <- sim[seq(1, (var *
                  rep), var) + (j - 1)]
            }
        }
    } else {
        stop("Family not known")
    }




    # Prepare data frame
    data <- data.frame(1:(rep * clus), paste0(rep("C", rep * clus),
        rep(1:clus, each = rep)))
    data <- data.frame(1:(rep * clus), paste0(rep("C", rep * clus),
        rep(1:clus, each = rep)), res)
    names(data)[1:2] <- c("Observation", "Cluster")

    return(data)
}
