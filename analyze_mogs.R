library(Hmisc)
library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringr)
library(ggthemes)

df_second_fitting_1 <- readRDS("df_second_fitting_1.rds")
df_second_fitting_2 <- readRDS("df_second_fitting_2.rds")
mogs <- readRDS("mogs.rds")

get.prop.skip <- function(position, target.mean, skip.mean,
                          sigma.target, sigma.skip,
                          lambda.target, lambda.skip) {
  ## gets proportion of intended skips given a particular landing site
  prob.targ <- pnorm(position+.5, target.mean, sigma.target) -
    pnorm(position-.5, target.mean, sigma.target)
  prob.skip <- pnorm(position+.5, skip.mean, sigma.skip) -
    pnorm(position-.5, skip.mean, sigma.skip)
  prop.skip <- (lambda.skip * prob.skip) /
    (lambda.skip * prob.skip + lambda.target * prob.targ)
  return(prop.skip)
}

mogs <- mogs %>%
  filter(sigma.skip >= 0.2, sigma.target >= 0.2, lambda.skip >= 0) %>%
  mutate(prop.skip.1 = get.prop.skip(1, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.2 = get.prop.skip(2, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.3 = get.prop.skip(3, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.4 = get.prop.skip(4, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.5 = get.prop.skip(5, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.6 = get.prop.skip(6, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.7 = get.prop.skip(7, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         diff14 = prop.skip.4 - prop.skip.1,
         diff25 = prop.skip.5 - prop.skip.2,
         diff36 = prop.skip.6 - prop.skip.3,
         diff47 = prop.skip.7 - prop.skip.4,
         diffmean = (diff14 + diff25 + diff36 + diff47) / 4)

mogs.1 <- mogs %>%
  filter(launch == 1)
mogs.2 <- mogs %>%
  filter(launch == 2)

df_plot_1 <- df_second_fitting_1 %>%
  filter(lpostarget <= 12)
df_plot_2 <- df_second_fitting_2 %>%
  filter(lpostarget <= 12)

get_max <- function(df, col_name) {
  d <- df[which.max(df[[col_name]]), ]
  d$maximizes <- col_name
  return(d)
}

max_mogs.1 <- bind_rows(lapply(c("diff14", "diff25", "diff36", "diff47", "loglik"),
                            function(x) get_max(mogs.1, x)))
max_mogs.2 <- bind_rows(lapply(c("diff14", "diff25", "diff36", "diff47", "loglik"),
                            function(x) get_max(mogs.2, x)))
max_mogs <- bind_rows(max_mogs.1, max_mogs.2) %>%
  mutate(label = str_c(maximizes, '_', launch))

df_second_fitting <- bind_rows(df_second_fitting_1, df_second_fitting_2) %>%
  mutate(launch = -launch)

print(max_mogs)

# my_colors <- scales::brewer_pal(type = "qual", pal = "Dark2")(8)
my_colors <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, 
                direction = 1)(4)
scales::show_col(my_colors)
tot_color <- my_colors[4]
skip_color <- my_colors[1]
targ_color <- my_colors[3]

## plot mixture-of-gaussians
plot_mog <- function(mog) {
  gaussian_target <- function(x) {
    dnorm(x, mean = mog$target.mean, sd = mog$sigma.target) * mog$lambda.target
  }
  gaussian_skip <- function(x) {
    dnorm(x, mean = mog$skip.mean, sd = mog$sigma.skip) * mog$lambda.skip
  }
  gaussian_total <- function(x) {
    gaussian_target(x) + gaussian_skip(x)
  }
  df <- df_second_fitting %>%
    filter(launch == mog$launch)
  b1 <- str_extract_all("ab_cruised_down", "[a-z_]")[[1]]
  b2 <- c(rep('', 3), as.character(seq(7)), rep('', 5))
  labels <- str_c(b1, b2, sep = '\n')
  p <- ggplot(df, aes(lpostarget)) +
    geom_histogram(breaks = seq(-2.5, 12.5), aes(y = ..density..)) +
    stat_function(fun = gaussian_target, color = targ_color) +
    stat_function(fun = gaussian_skip, color = skip_color) +
    stat_function(fun = gaussian_total, color = tot_color) +
    coord_cartesian(ylim = c(0, 0.31), xlim = c(-2.5, 12.5)) +
    scale_x_continuous(breaks = seq(-2, 12),
                       labels = labels) +
    labs(x = "Original destination",
         y = "Probability density")
#           panel.background = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())
  return(p)
}

max_mogs_as_list <- split(max_mogs, max_mogs$label)

mog_plots <- lapply(max_mogs_as_list, plot_mog)
# dummy <- mapply(function(plot, name) {
#   ggsave(str_c(name, "_mog.pdf"), plot, height = 8, width = 5.88)
# }, mog_plots, names(mog_plots))

## plot proportion of skipping saccades
plot_props <- function(mog) {
  b1 <- str_extract_all("cruised", "[a-z]")[[1]]
  b2 <- as.character(seq(7))
  labels <- str_c(b1, b2, sep = '\n')
  mogs <- mog %>%
    gather(lpos, prop_skip, prop.skip.1, prop.skip.2, prop.skip.3, prop.skip.4, prop.skip.5, prop.skip.6, prop.skip.7) %>%
    mutate(lpos = as.numeric(str_sub(lpos, -1, -1))) %>%
    mutate(population = "skip")

  p <- ggplot(mogs, aes(lpos, prop_skip, fill = skip_color)) + geom_bar(stat = "identity") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(7), 
                       labels = labels) +
    labs(x = "Original destination",
         y = "Proportion intended-skip") +
    theme(legend.position = "none")
  return(p)
}

prop_plots <- lapply(max_mogs_as_list, plot_props)
# dummy <- mapply(function(plot, name) {
#   ggsave(str_c(name, "_prop.pdf"), plot, height = 8, width = 5.88)
# }, prop_plots, names(prop_plots))

## plot log likelihood
plot_loglik <- function(mog) {
  p <- ggplot(mogs %>% filter(launch == mog$launch) %>%
                filter(loglik > quantile(loglik, probs = 0.5)), # top 20 % of models
              aes(launch, loglik)) + geom_violin(adjust = 0.7) +
    geom_hline(yintercept = mog$loglik, color = "red") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    labs(y = "Log-likelihood")
}

loglik_plots <- lapply(max_mogs_as_list, plot_loglik)

## plot everything together

pdf("mogs_1.pdf", height = 28/9, width = 7)
grid.arrange(
  mog_plots$loglik_1, prop_plots$loglik_1,
  ncol = 2, widths = c(4, 3) #, main = "title", left = "RT", sub = "LOG"
)
dev.off()

pdf("mogs_2.pdf", height = 28/9, width = 7)
grid.arrange(
  mog_plots$loglik_2, prop_plots$loglik_2,
  ncol = 2, widths = c(4, 3) #, main = "title", left = "RT", sub = "LOG"
)
dev.off()

pdf("mogs_1_SI.pdf", height = 15, width = 7.5)
grid.arrange(
  grobs = list(
    mog_plots$loglik_1, prop_plots$loglik_1, loglik_plots$loglik_1,
    mog_plots$diff14_1, prop_plots$diff14_1, loglik_plots$diff14_1,
    mog_plots$diff25_1, prop_plots$diff25_1, loglik_plots$diff25_1,
    mog_plots$diff36_1, prop_plots$diff36_1, loglik_plots$diff36_1,
    mog_plots$diff47_1, prop_plots$diff47_1, loglik_plots$diff47_1),
  ncol = 3, widths = c(1, 1, 0.5), sub = "Landing position on target" #, main = "title", left = "RT", sub = "LOG"
)
dev.off()

pdf("mogs_2_SI.pdf", height = 15, width = 7.5)
grid.arrange(
  grobs = list(
    mog_plots$loglik_2, prop_plots$loglik_2, loglik_plots$loglik_2,
    mog_plots$diff14_2, prop_plots$diff14_2, loglik_plots$diff14_2,
    mog_plots$diff25_2, prop_plots$diff25_2, loglik_plots$diff25_2,
    mog_plots$diff36_2, prop_plots$diff36_2, loglik_plots$diff36_2,
    mog_plots$diff47_2, prop_plots$diff47_2, loglik_plots$diff47_2),
  ncol = 3, widths = c(1, 1, 0.5), sub = "Landing position on target" #, main = "title", left = "RT", sub = "LOG"
)
dev.off()

## plot hypothetical data
mog_hyp <- data_frame(target.mean = 4.5, sigma.target = 1.5, lambda.target = 0.5, skip.mean = 7, sigma.skip = 1.7, lambda.skip = 0.5)

mog_hyp <- mog_hyp %>%
  filter(sigma.skip >= 0.2, sigma.target >= 0.2, lambda.skip >= 0) %>%
  mutate(prop.skip.1 = get.prop.skip(1, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.2 = get.prop.skip(2, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.3 = get.prop.skip(3, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.4 = get.prop.skip(4, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.5 = get.prop.skip(5, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.6 = get.prop.skip(6, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         prop.skip.7 = get.prop.skip(7, target.mean, skip.mean, sigma.target, sigma.skip, lambda.target, lambda.skip),
         diff14 = prop.skip.4 - prop.skip.1,
         diff25 = prop.skip.5 - prop.skip.2,
         diff36 = prop.skip.6 - prop.skip.3,
         diff47 = prop.skip.7 - prop.skip.4,
         diffmean = (diff14 + diff25 + diff36 + diff47) / 4)

prop_plot_hyp <- plot_props(mog_hyp)

dat <- data_frame(x = seq(-2.5, 12.5))

gaussian_target <- function(x) {
  dnorm(x, mean = mog_hyp$target.mean, sd = mog_hyp$sigma.target) * mog_hyp$lambda.target
}
gaussian_skip <- function(x) {
  dnorm(x, mean = mog_hyp$skip.mean, sd = mog_hyp$sigma.skip) * mog_hyp$lambda.skip
}
gaussian_total <- function(x) {
  gaussian_target(x) + gaussian_skip(x)
}

b1 <- str_extract_all("ab_cruised_down", "[a-z_]")[[1]]
b2 <- c(rep('', 3), as.character(seq(7)), rep('', 5))
labels <- str_c(b1, b2, sep = '\n')

mog_plot_hyp <- ggplot(dat, aes(x)) +
  stat_function(fun = gaussian_target, color = targ_color) +
  stat_function(fun = gaussian_skip, color = skip_color) +
  stat_function(fun = gaussian_total, color = tot_color) +
  coord_cartesian(ylim = c(0, 0.31), xlim = c(-2.5, 12.5)) +
  scale_x_continuous(breaks = seq(-2, 12),
                     labels = labels) +
  labs(x = "Original destination",
       y = "Probability density")

pdf("mogs_hyp.pdf", height = 28/9, width = 7)
grid.arrange(
  mog_plot_hyp, prop_plot_hyp,
  ncol = 2, widths = c(4, 3)
)
dev.off()
