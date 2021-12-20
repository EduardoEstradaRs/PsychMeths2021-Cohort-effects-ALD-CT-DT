## Fit models from
#  Estrada, Bunge, & Ferrer (2021).
#  Controlling for cohort effects in accelerated longitudinal
#  designs using continuous- and discrete-time dynamic models
#  Psychological Methods
#  https://doi.org/10.1037/met0000427

library(OpenMx)

## Load data ---------------------------------------

head(dW) # dW is the data set in wide format for SEM model
head(dL) # dL is the data set in list format for SSM in OpenMx (see Hunter, 2018)
#cases <- seq(1,nrow(dW)) # use either this line or the next one
cases <- seq(1,length(dL))

# Declare variables for SEM models (change this to match your data)
opmxL <- list()
opmxL$SEMvars$manifests <- paste0("oy", seq(0,14)) # observed variables
opmxL$SEMvars$latents <- paste0("y", seq(0,14))    # latent levels
opmxL$SEMvars$diff <- paste0("d", seq(1,14))       # latent changes
opmxL$measPaths <- mxPath(from=opmxL$SEMvars$latents, to=opmxL$SEMvars$manifests,
                          arrows=1, free=FALSE, values=1)
opmxL$measErr <- mxPath(from=c("y00","yA",opmxL$SEMvars$manifests),
                        arrows=2, free=TRUE, values=c(25,.7,rep(2, length(opmxL$SEMvars$manifests) )),
                        labels=c("y0v","yAv",rep("MerY", length(opmxL$SEMvars$manifests) )))


## Prepare invariant part of the OpenMx SSM models ------------------
#  Based on Hunter (2018)
opmxL$amat_ct <- mxMatrix(name = "A", "Full", 2, 2, free = c(T,F,F,F),
                          values = c(-.2,0,1,0),
                          dimnames = list( c("y0", "yA"), c("y0", "yA") ),
                          labels = c("b_y", NA,NA,NA),
                          lbound = c(-1, NA,NA,NA),
                          ubound = c(0, NA,NA,NA))

opmxL$bmat <- mxMatrix(name = "B", "Zero", 2, 1)

opmxL$cmat <- mxMatrix(name = "C", "Full", 1, 2, free = FALSE,
                       values = c(1,0), 
                       dimnames = list( c("y"), c("y0", "yA") ),
                       labels = c(NA, NA)  )

opmxL$dmat <- mxMatrix("Zero", 1, 1, name = "D")
opmxL$qmat <- mxMatrix("Zero", 2, 2, name = "Q")

opmxL$rmat <- mxMatrix("Diag", 1, 1, TRUE, 2,
                       name = "R", labels = "MerY")

# # Initial conditions predicted by covariates # #
opmxL$inicovar <- mxMatrix(name='iniCovars', nrow=2, ncol=1,
                           values=c(1, NA), labels=c(NA, 'data.coh'))

opmxL$inicoef <- mxMatrix(name='iniCoef', nrow=2, ncol=2,
                          values=c(10, 6, 0,0),
                          free = TRUE,
                          #free = c(T,T,F,T),
                          labels=c('y0mn', 'yAmn',
                                   'coh_y0', 'coh_yA'),
                          lbound = c(-100,-100,-100,-100),
                          ubound = c(100,100,100,100))

opmxL$xmat_reg <- with(opmxL, mxAlgebra(name='x0', iniCoef %*% iniCovars) )
# # # # # # # # # # # # # # # # # # # # #

opmxL$xmat <- mxMatrix(name = "x0", "Full", 2, 1, free = TRUE,
                       values = c(12, 7),
                       labels = c("y0mn", "yAmn"))

opmxL$pmat <- mxMatrix(name = "P0", "Symm", 2, 2, TRUE,
                       values = c(25, 3, .7),
                       labels = c("y0v", "y0Acv", "yAv"),
                       lbound = c(0, NA, 0))

# Next one must be a column vector with the same number of rows as the B and D matrices have columns
opmxL$umat <- mxMatrix("Zero", 1, 1, name = "u")

opmxL$tmat <- mxMatrix('Full', 1, 1, name='time', labels='data.age')


opmxL$modL_ct <- with(opmxL, list(amat_ct, bmat, cmat, dmat,
                                  qmat, rmat, xmat, pmat,
                                  umat, tmat))

opmxL$modL_ct_reg <- with(opmxL, list(amat_ct, bmat, cmat, dmat,
                                      qmat, rmat, xmat_reg, pmat,
                                      umat, tmat, inicoef, inicovar))

opmxL$expODE <- mxExpectationStateSpaceContinuousTime(A = "A", B = "B",
                                                      C = "C", D = "D",
                                                      Q = "Q", R = "R",
                                                      x0 = "x0", P0 = "P0",
                                                      u = "u", t = "time")

## Create SSM CT multisubject model
genMxIndModels_ct <- function(x, dwork, modNames_ct, coh=FALSE) {
  DataSetForSubjectK <- dwork[[x]]
  if(coh) {ctModel <- opmxL$modL_ct_reg} else {ctModel <- opmxL$modL_ct}
  indivmodels <- mxModel(name = modNames_ct[x],
                         ctModel,
                         opmxL$expODE,
                         mxFitFunctionML(),
                         mxData(DataSetForSubjectK, type ='raw')  )  }

  
##  Estimate SSM-CT without cohort effects --------------
modNames_ct <- paste0("i", cases, "ODE")
indivmodels_ct <- lapply(cases, genMxIndModels_ct, dL, modNames_ct)
multiSubjODE <- mxModel(name = "MultiODE", indivmodels_ct,
                        mxFitFunctionMultigroup(modNames_ct))

multiSubjODERun <- mxRun(multiSubjODE)
multiSubjODEsumm <- summary(multiSubjODERun)

## Estimate SSM-CT with cohort effects -----------
modNames_ct_coh <- paste0("i", cases, "ODEcoh")
indivmodels_ct <- lapply(cases, genMxIndModels_ct, dL, modNames_ct_coh, coh=TRUE)
multiSubjODEcoh <- mxModel(name = "MultiODEcoh", indivmodels_ct,
                           mxFitFunctionMultigroup(modNames_ct_coh))

multiSubjODEcohRun <- mxRun(multiSubjODEcoh)
multiSubjODEcohSumm <- summary(multiSubjODEcohRun)
      

## Estimate SEM without cohort effects -------------
LCSMx <- mxModel(name = "mxLCS_SEM", mxData(observed = dW, type="raw"),
                 type="RAM",
                 manifestVars = opmxL$SEMvars$manifests,
                 latentVars = c("y00", "yA", opmxL$SEMvars$latents, opmxL$SEMvars$diff),
                 opmxL$measPaths,
                 mxPath(from=opmxL$SEMvars$latents[1:14], to=opmxL$SEMvars$latents[2:15],
                        arrows=1, free=FALSE, values=1),
                 mxPath(from=opmxL$SEMvars$diff, to=opmxL$SEMvars$latents[2:15],
                        arrows=1, free=FALSE, values=1),
                 mxPath(from=opmxL$SEMvars$latents[1:14], to=opmxL$SEMvars$diff,
                        arrows=1, free=TRUE, values=-.2, labels="b_y"),
                 mxPath(from="y00", to="y0", arrows=1, free=FALSE, values=1),
                 mxPath(from="yA", to=opmxL$SEMvars$diff, arrows=1, free=FALSE, values=1),
                 mxPath(from="one", to=c("y00","yA"), free=TRUE, values=c(12,7),
                        labels=c("y0mn","yAmn")),
                 opmxL$measErr,
                 mxPath(from="y00", to="yA", arrows=2, free=TRUE, values=3, labels="y0Acv") )

SEMrun <- mxRun(LCSMx)  
SEMsumm <- summary(SEMrun)

## Estimate SEM with cohort effects -------------
LCSMxCoh <- mxModel(name = "mxLCS_SEMcoh", mxData(observed = dW, type="raw"),
                    type="RAM",
                    manifestVars = opmxL$SEMvars$manifests,
                    latentVars = c("y00", "yA", opmxL$SEMvars$latents, opmxL$SEMvars$diff, "cohort"),
                    
                    opmxL$measPaths,
                    opmxL$measErr,
                    
                    # Model dynamics
                    mxPath(from=opmxL$SEMvars$latents[1:14], to=opmxL$SEMvars$latents[2:15],
                           arrows=1, free=FALSE, values=1),
                    mxPath(from=opmxL$SEMvars$diff, to=opmxL$SEMvars$latents[2:15],
                           arrows=1, free=FALSE, values=1),
                    mxPath(from=opmxL$SEMvars$latents[1:14], to=opmxL$SEMvars$diff,
                           arrows=1, free=TRUE, values=-.2, labels="b_y"),
                    
                    # Latent initial state and additive component
                    mxPath(from="y00", to="y0", arrows=1, free=FALSE, values=1),
                    mxPath(from="yA", to=opmxL$SEMvars$diff, arrows=1, free=FALSE, values=1),
                    mxPath(from="one", to=c("y00","yA"), free=TRUE, values=c(12,7),
                           labels=c("y0mn","yAmn")),
                    mxPath(from="y00", to="yA", arrows=2, free=TRUE, values=3, labels="y0Acv"),
                    
                    # Definition variable / covariate
                    mxPath(from="one", to="cohort", free=FALSE, values=1,
                           labels="data.coh" ),
                    mxPath(from="cohort", to=c("y00","yA"), arrows=1,
                           free=TRUE, values=c(1,1), labels=c("coh_y0","coh_yA")  )
                    )

SEMCohRun <- mxRun(LCSMxCoh, intervals=TRUE)  
SEMCohRunSumm <- summary(SEMCohRun, verbose=T)


