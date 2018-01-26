library(doParallel)
library(mixtools)
library(plyr)
library(Hmisc)
library(dplyr)

source("mynormalmixEM.R")

registerDoParallel()

df.second <- readRDS("df.second.rds")

## mixture fitting
get.best.mog <- function(values, mean.constr, iterations) {
  maxrestarts <- 1000 # increase this if you get an error "too many tries!"
  maxit <- 10000 # increase this if you get an error "WARNING! NOT CONVERGENT!"
  failures <- rep(T, times=iterations) # pretend all models have failed
  models <- rep(list(NULL), times=iterations)
  first.run <- T
  while (any(failures)) {
    num.failures <- sum(failures)
    new.models <- replicate(num.failures,
                            normalmixEM(values,
                                        mean.constr=mean.constr,
                                        maxit=maxit,
                                        maxrestarts=maxrestarts),
                            simplify=F)
    models[failures] <- new.models # replace failed models
    failures <- (sapply(models, class) != "mixEM") # returns a mixEM object only
    # if it converges
    if (all(failures)) {
      stop("all", iterations, "models failed to converge.",
           "something is wrong!\n")
    } else if (first.run & any(failures)) {
      cat(sum(failures), "of", num.failures,
          "models failed to converge. Re-running them until convergence.\n")
      first.run <- F
    }
    ## models <- models[!failures]
  }
  likelihoods <- lapply(models, function(x) x$loglik)
  return(models[[which.max(likelihoods)]])
}

grid.size <- 0.1 # should use 0.1
outlier.cutoff <- seq(12, 12) # above position 12, there are no positions with > 5 observations in either dataset. position 12 is always above 5
target.means.1 <- seq(3.9, 5, grid.size)
target.means.2 <- seq(3.6, 5, grid.size)
skip.mean.diff <- seq(1.5, 3, grid.size)

df.mixes.1 <- expand.grid(target.mean=target.means.1, skip.mean.diff=skip.mean.diff, outlier.cutoff=outlier.cutoff)
df.mixes.1$skip.mean <- with(df.mixes.1, target.mean + skip.mean.diff)
df.mixes.2 <- expand.grid(target.mean=target.means.2, skip.mean.diff=skip.mean.diff, outlier.cutoff=outlier.cutoff)
df.mixes.2$skip.mean <- with(df.mixes.2, target.mean + skip.mean.diff)

df.second.fitting.1 <- df.second %>%
  filter(!is.na(lpos) & lpos>0 & distance >= 0 & launch == -1)
df.second.fitting.2 <- df.second %>%
  filter(!is.na(lpos) & lpos>0 & distance >= 0 & launch == -2)

get.two.constr.mean.mog.df <- function(values, target.mean,
                                       skip.mean, outlier.cutoff) {
  set.seed(11)
  iterations <- 100
  outliers <- (values > outlier.cutoff)
  v <- values[!outliers]
  mog <- get.best.mog(v,
                      mean.constr=c(target.mean, skip.mean), iterations)
  data.frame(loglik=mog$loglik,
             sigma.target=mog$sigma[1], sigma.skip=mog$sigma[2],
             lambda.target=mog$lambda[1], lambda.skip=mog$lambda[2])
}

mogs.1 <- ddply(df.mixes.1, .(target.mean, skip.mean, outlier.cutoff),
                function(df) {
                  get.two.constr.mean.mog.df(df.second.fitting.1$lpostarget,
                                             df$target.mean,
                                             df$skip.mean, df$outlier.cutoff)
                },
                .parallel=T)

mogs.2 <- ddply(df.mixes.2, .(target.mean, skip.mean, outlier.cutoff),
                function(df) {
                  get.two.constr.mean.mog.df(df.second.fitting.2$lpostarget,
                                             df$target.mean,
                                             df$skip.mean, df$outlier.cutoff)
                },
                .parallel=T)

mogs.1 <- mogs.1 %>%
  mutate(launch = 1)
mogs.2 <- mogs.2 %>%
  mutate(launch = 2)
mogs <- bind_rows(mogs.1, mogs.2)

saveRDS(mogs, "mogs.rds")
saveRDS(df.second.fitting.1, "df_second_fitting_1.rds")
saveRDS(df.second.fitting.2, "df_second_fitting_2.rds")
