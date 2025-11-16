#### NICOMPR Illustrations ####
#   by:         Meyer, MJ, Graye, A, & Sellers, KF
#   modified:   06/27/23
#

#### load crab data ####
# library(rsq) # data available on CRAN in rsq package
# data("hcrabs")

crabs       <- read.table('./Illustrations/crab.txt', header = FALSE)
names(crabs)  <- c('id', 'color', 'spine', 'weight', 'width', 'satellites')

#### load textile data ####
textile     <- read.table('./Illustrations/TextileFabrics.txt', header = TRUE)

#### hard-code word length data ####
hungarian   <- c(rep(1,1421), rep(2,12333), rep(3,20711), rep(4,15590), 
                 rep(5,5544), rep(6,1510), rep(7,289), rep(8,60), rep(9,1))

slovak      <- c(rep(1,7), rep(2,33), rep(3,49), rep(4,22), rep(5,6))

#### load code ####
source('Code/nicompr.R')
options(mc.cores = parallel::detectCores())

#### set up priors ####
x             <- c(2, 0)
priors_c1     <- list(a = 1, b = 1, c = 1)        # add one
priors_c2     <- list(a = sum(x), b = sum(log(factorial(x))), c = length(x)) # add two subjects, one with x_i = large, x_i = 0
priors_c01    <- list(a = 0.1, b = 0.1, c = 0.1)
priors_c001   <- list(a = 0.01, b = 0.01, c = 0.01)



#### crab satellites ####
X         <- crabs$satellites

##### conjugate models #####
crab_c1     <- suppressWarnings(cmp(X, priors = priors_c1,
                                    type = 'conjugate'))
crab_c2     <- suppressWarnings(cmp(X, priors = priors_c2,
                                    type = 'conjugate'))
crab_c01    <- suppressWarnings(cmp(X, priors = priors_c01,
                                    type = 'conjugate'))
crab_c001   <- suppressWarnings(cmp(X, priors = priors_c001,
                                    type = 'conjugate'))

##### flat model #####
crab_f      <- suppressWarnings(cmp(X, priors = list(a = 0, b = 0, c = 0),
                                    type = 'flat'))

##### Jeffreys' model #####
crab_j      <- suppressWarnings(cmp(X, type = 'Jeffreys',
                                    control = list(adapt_delta = 0.999,
                                                   max_treedepth = 15)))

summary(crab_c1)
summary(crab_c2)
summary(crab_c01)
summary(crab_c001)
summary(crab_f)
summary(crab_j)



#### textiles ####
Y         <- textile$NoOfFaults

##### conjugate models #####
text_c1     <- suppressWarnings(cmp(Y, priors = priors_c1,
                                    type = 'conjugate'))
text_c2     <- suppressWarnings(cmp(Y, priors = priors_c2,
                                    type = 'conjugate'))
text_c01    <- suppressWarnings(cmp(Y, priors = priors_c01,
                                    type = 'conjugate'))
text_c001   <- suppressWarnings(cmp(Y, priors = priors_c001,
                                    type = 'conjugate'))

##### flat model #####
text_f      <- suppressWarnings(cmp(Y, priors = list(a = 0, b = 0, c = 0),
                                    type = 'flat'))

##### Jeffreys' model #####
text_j      <- suppressWarnings(cmp(Y, type = 'Jeffreys',
                                    control = list(adapt_delta = 0.999,
                                                   max_treedepth = 15)))

summary(text_c1)
summary(text_c2)
summary(text_c01)
summary(text_c001)
summary(text_f)
summary(text_j)



#### Hungarian word length ####
H           <- hungarian

##### conjugate models #####
hun_c1    <- suppressWarnings(cmp(H, priors = priors_c1,
                                  type = 'conjugate'))
hun_c2    <- suppressWarnings(cmp(H, priors = priors_c2,
                                  type = 'conjugate'))
hun_c01   <- suppressWarnings(cmp(H, priors = priors_c01,
                                  type = 'conjugate'))
hun_c001  <- suppressWarnings(cmp(H, priors = priors_c001,
                                  type = 'conjugate'))

##### flat model #####
hun_f     <- suppressWarnings(cmp(H, priors = list(a = 0, b = 0, c = 0),
                                  type = 'flat'))

##### Jeffreys' model #####
hun_j     <- suppressWarnings(cmp(H, type = 'Jeffreys',
                                  control = list(adapt_delta = 0.999,
                                                 stepsize = 0.01,
                                                 max_treedepth = 15)))
summary(hun_c1)
summary(hun_c2)
summary(hun_c01)
summary(hun_c001)
summary(hun_f)
summary(hun_j)



#### Slovak word length ####
S           <- slovak

##### conjugate models #####
slo_c1    <- suppressWarnings(cmp(S, priors = priors_c1,
                                  type = 'conjugate'))
slo_c2    <- suppressWarnings(cmp(S, priors = priors_c2,
                                  type = 'conjugate'))
slo_c01   <- suppressWarnings(cmp(S, priors = priors_c01,
                                  type = 'conjugate'))
slo_c001  <- suppressWarnings(cmp(S, priors = priors_c001,
                                  type = 'conjugate'))

##### flat model #####
slo_f     <- suppressWarnings(cmp(S, priors = list(a = 0, b = 0, c = 0),
                                  type = 'flat'))

##### Jeffreys' model #####
slo_j     <- suppressWarnings(cmp(S, type = 'Jeffreys',
                                  control = list(adapt_delta = 0.999,
                                                 stepsize = 0.01,
                                                 max_treedepth = 15)))

summary(slo_c1)
summary(slo_c2)
summary(slo_c01)
summary(slo_c001)
summary(slo_f)
summary(slo_j)

