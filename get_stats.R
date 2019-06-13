library(lme4)
library(MASS)
library(dplyr)
library(tidyr)
library(stringr)
library(parallel)

if (!require('lmer.slimr')) {
  devtools::install_github('kbicknell/lmer.slimr')
  library('lmer.slimr')
}

options(width = 160L)

df_lmer <- readRDS("df.lmer.rds")


full_dataset <- TRUE

if (full_dataset) {
  launch_sites <- c("launch3.pre2.1", "launch3.pre3.2", "launch3.pre4.3",
                    "launch3.pre5.4", "launch3.pre6.5", "launch3.pre7.6",
                    "launch3.pre8.7", "launch3.pre9.8", "launch3.pre10.9", "launch3.pre11.10")
} else {
  launch_sites <- "launch3.pre2.1"
}

shifts <- c(None = 0, Right = 3, Left = -3)
df_lmer <- df_lmer %>%
  mutate(shift = shifts[cond],
         population = factor(lpos3c + shift))

lmer_control <- lmerControl(optCtrl = list(maxfun = 100000L))
glmer_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L))

add_arb_contrasts <- function(df,
                              interpretations, # a matrix with named rows
                              old_vars, # a character vector of strings
                              other_vars_to_show, # a character vector of strings
                              other_vars_dont_show, # a character vector of strings
                              formula # just the fixed effects
) {
  df_contrasts <- df %>%
    group_by_(.dots = old_vars) %>%
    summarize() %>%
    ungroup()
  if (nrow(df_contrasts) != ncol(interpretations)) {
    print(df_contrasts)
    stop("Number of columns in interpretations must match unique combinations",
         "of old_vars (see just-printed dataframe)")
  }
  contrast_matrix <- ginv(interpretations)
  colnames(contrast_matrix) <- rownames(interpretations)
  rownames(contrast_matrix) <- do.call(paste, df_contrasts[old_vars])
  contrast_id <- do.call(paste, df[old_vars])
  for (contrast_name in rownames(interpretations)) {
    df[[contrast_name]] <- contrast_matrix[match(contrast_id, rownames(contrast_matrix)), contrast_name]
  }

  ## ensure that this contrast matrix meets the assumptions
  if (rankMatrix(cbind(1, contrast_matrix)) != nrow(contrast_matrix)) {
    warning("Contrast matrix is not of full rank.",
            "Only use this if you know what you're doing!")
  }
  unit_vector <- rep(sqrt(1 / nrow(contrast_matrix)), times = nrow(contrast_matrix))
  dim(unit_vector) <- c(nrow(contrast_matrix), 1)
  if (!(isTRUE(all.equal(unit_vector, Null(contrast_matrix))) ||
          isTRUE(all.equal(-unit_vector, Null(contrast_matrix))))) {
    warning("Linear restriction on the untransformed \\alpha_j values",
            "is not the usual assumption that \\sum_j \\alpha_j == 0.",
            "Instead, the basis vector of the null space is",
            as.character(Null(contrast_matrix)))
  }
  ## test whether this works with other factors involved
  df_test <- df %>%
    group_by_(.dots = c(old_vars, other_vars_to_show, other_vars_dont_show, rownames(interpretations))) %>%
    summarize(V1 = 1)
  for (old_var in rev(old_vars)) {
    df_test <- df_test[order(df_test[[old_var]]), ]
  }
  f <- formula(formula)
  m <- fractions(ginv(model.matrix(f, df_test)))
  
  rownames(m) <- names(coef(lm(paste0("V1", formula), df_test)))
  colnames(m) <- do.call(paste, c(df_test[c(old_vars, other_vars_to_show)], sep="/"))
  cat("Verify that this is the correct interpretation of variables.\nColumn headings show: ")
  cat(paste(c(old_vars, other_vars_to_show), collapse="/"), "\n")
  print(m)
  return(df)
}

#' gets gaze duration and refixation rate lme4 model
#' 
#' @param header a string to print at the top
#' @param expts a character vector of the experiments to include, e.g., c('yd', 'ydr')
#' @param interpretations a matrix with row names giving interpretation of each predictor
#' @param ranef_vars a character vector of variables to include in random effects
#' @param old_vars a character vector of variables to pass on to add_arb_contrasts
#' @param other_vars_to_show a character vector that will be passed on to add_arb_contrasts and included in the regression
#' @param other_vars_dont_show a character vector of variables to pass on to add_arb_contrasts
#' 
get_models <- function(header, df, interpretations, ranef_vars, old_vars,
                       other_vars_to_show, other_vars_dont_show) {
  cat(sprintf("\n\n:: %s ::\n", header))
  df <- df %>%
    mutate(lpos3c = factor(lpos3c),
           launch3.pre = factor(launch3.pre))
  launch3_pre_cols <- sdif.split(df$launch3.pre, "launch3.pre")
  df <- cbind(df, launch3_pre_cols)
  
  fixed_effects <- str_c(c(other_vars_to_show, row.names(interpretations)),
                         collapse = ' + ')
  fixed_effects <- paste0('~ ', fixed_effects)
  ranef_vars_str <- str_c(ranef_vars, collapse = ' + ')
  random_effects <- sprintf("(%s | subj) + (%s | item)", ranef_vars_str, ranef_vars_str)
  df <- add_arb_contrasts(df, interpretations, old_vars,
                          other_vars_to_show, other_vars_dont_show,
                          fixed_effects)
  
  rhs <- str_c(fixed_effects, ' + ', random_effects)
  gzd_formula <- as.formula(str_c("gzd3c", rhs))
  refix_formula <- as.formula(str_c("refix3c", rhs))
  environment(gzd_formula) <- environment() # apparently required by profile and sometimes drop1
  environment(refix_formula) <- environment()
  p_gzd <-
    mcparallel(lmer.slimr(gzd_formula, df,
                          simplest.only = F,
                          return.only.converged = T,
                          control = lmer_control))
  p_refix <-
    mcparallel(glmer.slimr(refix_formula, df, family = "binomial",
                           simplest.only = F,
                           return.only.converged = T,
                           control = glmer_control))
  m_gzd <- mccollect(p_gzd)[[1]]
  m_refix <- mccollect(p_refix)[[1]]
  return(list(gzd = m_gzd, refix = m_refix))
}

#' like drop1, but does a single LRT test removing a number of predictors, and returns the anova output
drop_many <- function(m, vars) {
  ## get the original dataframe from its environment...
  env <- environment(formula(m))
  vars_string <- str_c('- ', vars, collapse = ' ')
  update_string <- str_c('. ~ . ', vars_string)
  nfit <- update(m, as.formula(update_string), evaluate = FALSE)
  m_simpler <- eval(nfit, envir = env)
  anova_output <- anova(m, m_simpler)
  return(anova_output)
}

#' get the total LRT for removing all interaction terms in main effect models
get_interaction_lrts <- function(m) {
  f <- as.character(formula(m))[3]
  if (str_detect(f, 'i3')) {
    vars <- c("i1", "i2", "i3")
  } else {
    vars <- c("i1", "i2")
  }
  return(drop_many(m, vars))
}

#' gets the results of likelihood ratio tests in parallel and rbinds results
get_lrts <- function(m, vars_to_test) {
  lrts <- mclapply(vars_to_test,
                   function(v) {
                     drop1(m, scope = v, test = "Chisq")
                   })
  return(do.call(rbind, lrts))
}

#' fit two mixed-effects models to yosdos style dataset: one with main effect
#' and one without. analyze both with Wald and LRT tests
analyze_four_ways <- function(label, df, ints_me, ints_no_me, other_vars_to_show, no7 = FALSE) {

  me_random_effects <- c("odiff", "i1", "i2", "i3")
  if (no7) {
    me_random_effects <- me_random_effects[1:3]
  }
  
  models_me <- get_models(
    str_c(label, ": with main effect"),
    df,
    ints_me,
    me_random_effects,
    c("lpos3c", "population"),
    other_vars_to_show,
    c()
  )

  vars_to_test <- me_random_effects
  lrts_me <- mclapply(models_me, function(m) {
      get_lrts(m, vars_to_test)
  })
  lrts_int_me <- mclapply(models_me, function(m) {
    get_interaction_lrts(m)
  })
  
  print(summary(models_me[["gzd"]]))
  print(lrts_me[["gzd"]])
  print(lrts_int_me[["gzd"]])
  
  print(summary(models_me[["refix"]]))
  print(lrts_me[["refix"]])
  print(lrts_int_me[["refix"]])
  
  no_me_random_effects <- c("c4vs1", "c5vs2", "c6vs3", "c7vs4")
  if (no7) {
    no_me_random_effects <- no_me_random_effects[1:3]
  }
  models_no_me <- get_models(
    str_c(label, ": no main effect"),
    df,
    ints_no_me,
    no_me_random_effects,
    c("lpos3c", "population"),
    other_vars_to_show,
    c()
  )
  
  vars_to_test <- no_me_random_effects
  lrts_no_me <- mclapply(models_no_me, function(m) {
      get_lrts(m, vars_to_test)
  })
  
  print(summary(models_no_me[["gzd"]]))
  print(lrts_no_me[["gzd"]])
  
  print(summary(models_no_me[["refix"]]))
  print(lrts_no_me[["refix"]])

  
  return(list(me = list(models = models_me, lrts = lrts_me, lrts_int = lrts_int_me),
              no_me = list(models = models_no_me, lrts = lrts_no_me)))
}

analyze_yd_style_no7 <- function(label, df, other_vars_to_show = launch_sites) {
  df <- df %>%
    filter(lpos3c != 7, population != 7)
  df <- df %>%
    mutate(population = factor(population))
  
  ## intermediates to specify the main effect models easily
  c4vs1 <- c(-1/2, 1/2, 0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0)
  c5vs2 <- c(0, 0, -1/2, 1/2, 0, 0, 0, 0, -1/2, 1/2, 0, 0)
  c6vs3 <- c(0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0, -1/2, 1/2)
  c4vs1rml <- c(1/2, -1/2, 0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0)
  c5vs2rml <- c(0, 0, 1/2, -1/2, 0, 0, 0, 0, -1/2, 1/2, 0, 0)
  c6vs3rml <- c(0, 0, 0, 0, 1/2, -1/2, 0, 0, 0, 0, -1/2, 1/2)
  lpos2.1 <- c(-1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0, 0, 0)
  lpos3.2 <- c(0, 0, -1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0)
  lpos4.3 <- c(0, 0, 0, 0, -1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0)
  lpos5.4 <- c(0, 0, 0, 0, 0, 0, -1/2, -1/2, 1/2, 1/2, 0, 0)
  lpos6.5 <- c(0, 0, 0, 0, 0, 0, 0, 0, -1/2, -1/2, 1/2, 1/2)
  
  ints_me <- rbind(
    odiff = (c4vs1 + c5vs2 + c6vs3) / 3,
    i1 = c6vs3 - c4vs1,
    i2 = c6vs3 - c5vs2,
    c4vs1rml = c4vs1rml,
    c5vs2rml = c5vs2rml,
    c6vs3rml = c6vs3rml,
    lpos2.1 = lpos2.1,
    lpos3.2 = lpos3.2,
    lpos4.3 = lpos4.3,
    lpos5.4 = lpos5.4,
    lpos6.5 = lpos6.5
  )
  
  ints_no_me <- rbind(
    c4vs1 = c4vs1,
    c5vs2 = c5vs2,
    c6vs3 = c6vs3,
    c4vs1rml = c4vs1rml,
    c5vs2rml = c5vs2rml,
    c6vs3rml = c6vs3rml,
    lpos2.1 = lpos2.1,
    lpos3.2 = lpos3.2,
    lpos4.3 = lpos4.3,
    lpos5.4 = lpos5.4,
    lpos6.5 = lpos6.5
  )
  
  return(analyze_four_ways(label, df, ints_me, ints_no_me, other_vars_to_show, no7 = TRUE))
}

analyze_yd_style <- function(label, df, other_vars_to_show = launch_sites) {
  ## intermediates to specify the main effect models easily
  c4vs1 <- c(-1/2, 1/2, 0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0, 0, 0, 0)
  c5vs2 <- c(0, 0, -1/2, 1/2, 0, 0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0)
  c6vs3 <- c(0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0, 0, -1/2, 1/2, 0, 0)
  c7vs4 <- c(0, 0, 0, 0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0, -1/2, 1/2)
  c4vs1rml <- c(1/2, -1/2, 0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0, 0, 0, 0)
  c5vs2rml <- c(0, 0, 1/2, -1/2, 0, 0, 0, 0, 0, -1/2, 1/2, 0, 0, 0, 0)
  c6vs3rml <- c(0, 0, 0, 0, 1/2, -1/2, 0, 0, 0, 0, 0, -1/2, 1/2, 0, 0)
  c7vs4rml <- c(0, 0, 0, 0, 0, 0, 0, 1/2, -1/2, 0, 0, 0, 0, -1/2, 1/2)
  lpos2.1 <- c(-1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  lpos3.2 <- c(0, 0, -1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  lpos4.3 <- c(0, 0, 0, 0, -1/2, -1/2, 1/3, 1/3, 1/3, 0, 0, 0, 0, 0, 0)
  lpos5.4 <- c(0, 0, 0, 0, 0, 0, -1/3, -1/3, -1/3, 1/2, 1/2, 0, 0, 0, 0)
  lpos6.5 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1/2, -1/2, 1/2, 1/2, 0, 0)
  lpos7.6 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/2, -1/2, 1/2, 1/2)
  
  ints_me <- rbind(
    odiff = (c4vs1 + c5vs2 + c6vs3 + c7vs4) / 4,
    i1 = c7vs4 - c4vs1,
    i2 = c7vs4 - c5vs2,
    i3 = c7vs4 - c6vs3,
    c4vs1rml = c4vs1rml,
    c5vs2rml = c5vs2rml,
    c6vs3rml = c6vs3rml,
    c7vs4rml = c7vs4rml,
    lpos2.1 = lpos2.1,
    lpos3.2 = lpos3.2,
    lpos4.3 = lpos4.3,
    lpos5.4 = lpos5.4,
    lpos6.5 = lpos6.5,
    lpos7.6 = lpos7.6
  )
  
  ints_no_me <- rbind(
    c4vs1 = c4vs1,
    c5vs2 = c5vs2,
    c6vs3 = c6vs3,
    c7vs4 = c7vs4,
    c4vs1rml = c4vs1rml,
    c5vs2rml = c5vs2rml,
    c6vs3rml = c6vs3rml,
    c7vs4rml = c7vs4rml,
    lpos2.1 = lpos2.1,
    lpos3.2 = lpos3.2,
    lpos4.3 = lpos4.3,
    lpos5.4 = lpos5.4,
    lpos6.5 = lpos6.5,
    lpos7.6 = lpos7.6
  )
  
  return(analyze_four_ways(label, df, ints_me, ints_no_me, other_vars_to_show))
}

analyze_yos_style_no7 <- function(label, df, other_vars_to_show = launch_sites) {
  df <- df %>%
    filter(lpos3c <= 3, population != 7)
  df <- df %>%
    mutate(population = factor(population))
  
  ## intermediates to specify the main effect models easily
  c4vs1 <- c(-1, 1, 0, 0, 0, 0)
  c5vs2 <- c(0, 0, -1, 1, 0, 0)
  c6vs3 <- c(0, 0, 0, 0, -1, 1)
  lpos2.1 <- c(-1/2, -1/2, 1/2, 1/2, 0, 0)
  lpos3.2 <- c(0, 0, -1/2, -1/2, 1/2, 1/2)
  
  ints_me <- rbind(
    odiff = (c4vs1 + c5vs2 + c6vs3) / 3,
    i1 = c6vs3 - c4vs1,
    i2 = c6vs3 - c5vs2,
    lpos2.1 = lpos2.1,
    lpos3.2 = lpos3.2
  )
  
  ints_no_me <- rbind(
    c4vs1 = c4vs1,
    c5vs2 = c5vs2,
    c6vs3 = c6vs3,
    lpos2.1 = lpos2.1,
    lpos3.2 = lpos3.2
  )
  
  return(analyze_four_ways(label, df, ints_me, ints_no_me, other_vars_to_show, no7 = TRUE))
}

analyze_yos_style <- function(label, df, other_vars_to_show = launch_sites) {
  ## intermediates to specify the main effect models easily
  c4vs1 <- c(-1, 1, 0, 0, 0, 0, 0, 0)
  c5vs2 <- c(0, 0, -1, 1, 0, 0, 0, 0)
  c6vs3 <- c(0, 0, 0, 0, -1, 1, 0, 0)
  c7vs4 <- c(0, 0, 0, 0, 0, 0, -1, 1)
  lpos2.1 <- c(-1/2, -1/2, 1/2, 1/2, 0, 0, 0, 0)
  lpos3.2 <- c(0, 0, -1/2, -1/2, 1/2, 1/2, 0, 0)
  lpos4.3 <- c(0, 0, 0, 0, -1/2, -1/2, 1/2, 1/2)
  
  ints_me <- rbind(
    odiff = (c4vs1 + c5vs2 + c6vs3 + c7vs4) / 4,
    i1 = c7vs4 - c4vs1,
    i2 = c7vs4 - c5vs2,
    i3 = c7vs4 - c6vs3,
    lpos2.1 = lpos2.1,
    lpos3.2 = lpos3.2,
    lpos4.3 = lpos4.3
  )
  
  ints_no_me <- rbind(
    c4vs1 = c4vs1,
    c5vs2 = c5vs2,
    c6vs3 = c6vs3,
    c7vs4 = c7vs4,
    lpos2.1 = lpos2.1,
    lpos3.2 = lpos3.2,
    lpos4.3 = lpos4.3
  )
  
  return(analyze_four_ways(label, df, ints_me, ints_no_me, other_vars_to_show))
}

#' @param type either 'yos' or 'yd' style
get_results <- function(type, label, df, other_vars_to_show = launch_sites) {
  results <- list()
  
  if (type == 'yd') {
    analyze_fn <- analyze_yd_style
    analyze_fn_no7 <- analyze_yd_style_no7
  } else if (type == 'yos') {
    analyze_fn <- analyze_yos_style
    analyze_fn_no7 <- analyze_yos_style_no7
  } else {
    stop("invalid type", type)
  }
  
  p_full <- (analyze_fn(label, df, other_vars_to_show))
  p_no7 <- (analyze_fn_no7(paste0(label, ": no 7"), df, other_vars_to_show))
  
  ## results$full <- mccollect(p_full)[[1]]
  ## results$no7 <- mccollect(p_no7)[[1]]
  results$full <- p_full
  results$no7 <- p_no7
  return(results)
}

results_full <- get_results(
  'yd',
  "Full dataset",
  df_lmer)

results_ydydr <- get_results(
  'yd',
  "YosDos and YosDosRep",
  df_lmer %>% filter(expt %in% c("yd", "ydr")))

results_ydr <- get_results(
  'yd',
  "YosDosRep",
  df_lmer %>% filter(expt == "ydr"))

results_yd <- get_results(
  'yd',
  "YosDos",
  df_lmer %>% filter(expt == "yd"))

# this analysis seems to be broken since launch site doesn't vary for lpos3c == 1, population == 1
results_yos_broken <- get_results(
  'yos',
  "Yos: with launch site (broken)",
  df_lmer %>% filter(expt == "yos", lpos3c %in% c(1, 2, 3, 4)))

results_yos_no_launch <- get_results(
  'yos',
  "Yos: without launch site",
  df_lmer %>% filter(expt == "yos", lpos3c %in% c(1, 2, 3, 4)),
  c()) # don't include launch site

save.image("stats.RData")
