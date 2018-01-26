# let's visualize CIs!
library(Hmisc)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)

load("stats.RData")

results_all <- list()
results_all$full <- results_full
results_all$ydydr <- results_ydydr
results_all$ydr <- results_ydr
results_all$yd <- results_yd
results_all$yos <- results_yos_no_launch
# results_all$yos_broken <- results_yos_broken


#' returns a data frame of CIs from a merMod
get_CIs <- function(m) {
  if (!isClass(m, "merMod")) { # sometimes no model here (when it didn't converge)
    return(data.frame())
  }
  df <- as.data.frame(confint(m, parm = "beta_", method = "Wald"))
  df$variable <- row.names(df)
  df$beta <- m@beta
  df <- df %>%
    select(lower = starts_with("2.5"),
           upper = starts_with("97.5"),
           variable, beta)
  return(df)
}

m <- results_all$full$full$me$models$refix

blah <- results_all$full

get_model_CIs <- function(x) {
  return(ldply(x, get_CIs, .id = "measure"))
}

get_CIs_from_models <- function(x) {
  return(get_model_CIs(x$models))
}

get_me_and_no_me <- function(x) {
  ldply(x, get_CIs_from_models, .id = "me")
}

get_7_and_no_7 <- function(x) {
  ldply(x, get_me_and_no_me, .id = "seven")
}

get_all_expts <- function(x) {
  ldply(x, get_7_and_no_7, .id = 'expt')
}

CIs <- get_all_expts(results_all)

ggplot(CIs %>% filter(variable %in% c('odiff', 'c4vs1', 'c5vs2', 'c6vs3', 'c7vs4'), measure == 'gzd'), aes(x = expt, y = beta)) + geom_pointrange(aes(ymin = lower, ymax = upper)) + coord_flip() + facet_grid(seven ~ variable)
ggsave("gzd_all_expts.pdf", height = 4, width = 10)

ggplot(CIs %>% filter(variable %in% c('odiff', 'c4vs1', 'c5vs2', 'c6vs3', 'c7vs4'), measure == 'refix'), aes(x = expt, y = beta)) + geom_pointrange(aes(ymin = lower, ymax = upper)) + facet_grid(seven ~ variable) + coord_flip() + ylim(-3, 2)
ggsave("refix_all_expts.pdf", height = 4, width = 10)

CIs_overall <- CIs %>%
  filter(variable %in% c("odiff", "c4vs1", "c5vs2", "c6vs3", "c7vs4"),
         expt == 'full', seven == 'full') %>%
  mutate(variable = factor(variable, levels = c("odiff", "c4vs1", "c5vs2", "c6vs3", "c7vs4")),
         beta = -beta,
         upper = -upper,
         lower = -lower)
levels(CIs_overall$variable) <- c("Overall\naverage", "1 vs. 4", "2 vs. 5", "3 vs. 6", "4 vs. 7")

ggplot(CIs_overall %>% filter(measure == 'gzd'), aes(x = variable, y = beta)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  theme(axis.title.x = element_blank()) +
    labs(y = "Coefficient estimate (ms)") +
    coord_cartesian(ylim=c(-20, 70))
ggsave("gzd_betas_overall.pdf", height = 3.5, width = 3.5)

ggplot(CIs_overall %>% filter(measure == 'refix'), aes(x = variable, y = beta)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  theme(axis.title.x = element_blank()) +
    labs(y = "Coefficient estimate (logits)") +
    coord_cartesian(ylim=c(-0.25, 2.5))
ggsave("refix_betas_overall.pdf", height = 3.5, width = 3.5)

