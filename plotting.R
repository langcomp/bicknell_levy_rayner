## PLOTTING
library(ggplot2)
library(ggthemes)
library(Hmisc)
library(dplyr)
library(tidyr)
library(stringr)

if (!require('ggintervals')) {
  devtools::install_github('langcomp/ggintervals')
  library('ggintervals')
}

df_lmer <- readRDS("df.lmer.rds")

df_plot <- df_lmer %>%
  group_by(subj, label, lpos3c, cond, expt, launch3.pre) %>%
  summarize(gzd3c = mean(gzd3c),
            refix3c = mean(refix3c),
            single3c = mean(single3c, na.rm = TRUE),
            ffix3c = mean(ffix3c)) %>%
  ungroup()

simplify_all_same_vector <- function(x) {
  if (all(x == x[1])) {
    return(x[1])
  } else {
    stop("x not all equal", x)
  }
}

## fancy new graph
comp_left_break <- "Further\nbehind"
comp_right_break <- "Further\nforward"
comp_left <- "Further behind"
comp_right <- "Further forward"
population_name <- "Original destination"
lpos_name <- "Actual landing site"

get_panel <- function(lpos3c, cond) {
  lpos3c <- simplify_all_same_vector(lpos3c)
  cond <- simplify_all_same_vector(cond)
  if (lpos3c %in% c(2, 5)) {
    return("2 vs. 5")
  } else if (lpos3c %in% c(3, 6)) {
    return("3 vs. 6")
  } else if (lpos3c == 1) {
    return("1 vs. 4")
  } else if (lpos3c == 7) {
    return("4 vs. 7")
  } else if (lpos3c == 4) {
    if (cond %in% c("None", "Left")) {
      return("1 vs. 4")
    } else if (cond == "Right") {
      return("4 vs. 7")
    }
  } else {
    stop("invalid lpos3c", lpos3c)
  }
}

df_fancyplot <- df_plot %>%
  group_by(lpos3c, cond) %>%
  mutate(panel = get_panel(lpos3c, cond)) %>%
  ungroup()

## duplicate lpos 4 no shift data for panel 4 vs 7 too
df_fancyplot_4n <- df_fancyplot %>%
  filter(lpos3c == 4, cond == "None") %>%
  mutate(panel = "4 vs. 7")

df_fancyplot <- rbind(df_fancyplot, df_fancyplot_4n) %>%
  mutate(panel_left = as.numeric(str_sub(panel, 1, 1)),
         comparison = ifelse(panel_left == lpos3c, comp_left_break, comp_right_break),
         population = ifelse(((comparison == comp_left_break) & (cond == "None")) |
                               ((comparison == comp_right_break) & (cond == "Left")),
                             comp_left, comp_right))

mean_cl_boot_bca_10k <- function(...) {mean_cl_boot_bca(..., replicates = 10000)}
## without launch site
gzd_nolaunch_plot <- ggplot(df_fancyplot, aes(as.factor(lpos3c), gzd3c, color=population, group=population)) +
  stat_summary(fun.data="mean_cl_boot_bca_10k", position=position_dodge(width=.3), fun.args = list(conf.int=.95)) +
  facet_grid(~ panel, scale = "free_x") +
  labs(x = lpos_name, y="Gaze duration (ms)", color = population_name,
       title = "Pair of original destinations being compared") +
    theme(plot.title = element_text(size=12)) +
    scale_color_colorblind()
ggsave("gzd_nolaunch_full_ext.pdf", gzd_nolaunch_plot, height=4, width=7)

refix_nolaunch_plot <- ggplot(df_fancyplot, aes(as.factor(lpos3c), refix3c, color=population, group=population)) +
  stat_summary(fun.data="mean_cl_boot_bca_10k", fun.args = list(conf.int=.95), position=position_dodge(width=.3)) +
  facet_grid(~ panel, scale = "free_x") +
  labs(x = lpos_name, y="Refixation probability", color = population_name,
       title = "Pair of original destinations being compared") +
  theme(plot.title = element_text(size=12)) +
    scale_color_colorblind()
ggsave("refix_nolaunch_full_ext.pdf", refix_nolaunch_plot, height=4, width=7)

sfd_nolaunch_plot <- ggplot(df_fancyplot, aes(as.factor(lpos3c), single3c, color=population, group=population)) +
  stat_summary(fun.data="mean_cl_boot_bca_10k", position=position_dodge(width=.3), fun.args = list(conf.int=.95)) +
  facet_grid(~ panel, scale = "free_x") +
  labs(x = lpos_name, y="Single fixation duration (ms)", color = population_name,
       title = "Pair of original destinations being compared") +
  theme(plot.title = element_text(size=12)) +
    scale_color_colorblind()
ggsave("sfd_nolaunch_full_ext.pdf", sfd_nolaunch_plot, height=4, width=7)

ffd_nolaunch_plot <- ggplot(df_fancyplot, aes(as.factor(lpos3c), ffix3c, color=population, group=population)) +
  stat_summary(fun.data="mean_cl_boot_bca_10k", position=position_dodge(width=.3), fun.args = list(conf.int=.95)) +
  facet_grid(~ panel, scale = "free_x") +
  labs(x = lpos_name, y="First fixation duration (ms)", color = population_name,
       title = "Pair of original destinations being compared") +
  theme(plot.title = element_text(size=12)) +
    scale_color_colorblind()
ggsave("ffd_nolaunch_full_ext.pdf", ffd_nolaunch_plot, height=4, width=7)

## with launch site (without bca bootstrap method, since bca fails in certain cases)
gzd_launch_plot <- ggplot(df_fancyplot, aes(as.factor(lpos3c), gzd3c, color=population, group=population)) +
  stat_summary(fun.data="mean_cl_boot",  fun.args = list(conf.int=.95), position=position_dodge(width=.3)) +
  facet_grid(launch3.pre ~ panel, scale = "free_x") +
  labs(x = lpos_name, y="Gaze duration (ms)", color = population_name) +
    scale_color_colorblind()
ggsave("gzd_launch_full.pdf", gzd_launch_plot, height=6, width=7)

refix_launch_plot <- ggplot(df_fancyplot, aes(as.factor(lpos3c), refix3c, color=population, group=population)) +
  stat_summary(fun.data="mean_cl_boot",  fun.args = list(conf.int=.95), position=position_dodge(width=.3)) +
  facet_grid(launch3.pre ~ panel, scale = "free_x") +
  labs(x = lpos_name, y="Refixation probability", color = population_name) +
    scale_color_colorblind()
ggsave("refix_launch_full.pdf", refix_launch_plot, height=6, width=7)


# ## if we want diffs, we have to do them between-subjects (because of data sparsity)
# df.fancy.diffs <- ddply(df.fancyplot, .(comparison, panel), function(df) {
#   test.gzd <- with(df, t.test(gzd3c ~ population))
#   ymin.gzd <- test.gzd$conf.int[1]
#   ymax.gzd <- test.gzd$conf.int[2]
#   y.gzd <- mean(c(ymin.gzd, ymax.gzd))
#   test.refix <- with(df, t.test(refix3c ~ population))
#   ymin.refix <- test.refix$conf.int[1]
#   ymax.refix <- test.refix$conf.int[2]
#   y.refix <- mean(c(ymin.refix, ymax.refix))
#   lpos3c <- df$lpos3c[1]
#   panel <- df$panel[1]
#   panel.left <- df$panel.left[1]
#   comparison <- df$comparison[1]
#   return(data.frame(ymin.gzd=ymin.gzd, ymax.gzd=ymax.gzd, y.gzd=y.gzd,
#                     ymin.refix=ymin.refix, ymax.refix=ymax.refix, y.refix=y.refix,
#                     lpos3c=lpos3c, panel=panel, panel.left=panel.left,
#                     comparison=comparison))
# })
# 
# ggplot(df.fancy.diffs, aes(x=comparison)) + geom_pointrange(aes(y=y.gzd, ymin=ymin.gzd, ymax=ymax.gzd)) + geom_hline(aes(y=0, color='dummy')) + facet_wrap(~panel)
# ggplot(df.fancy.diffs, aes(x=comparison)) + geom_pointrange(aes(y=y.refix, ymin=ymin.refix, ymax=ymax.refix)) + geom_hline(aes(y=0, color='dummy')) + facet_wrap(~panel)
# 
# ggplot(df.fancy.diffs, aes(x=lpos3c, color=comparison, group=comparison)) + geom_pointrange(aes(y=y.gzd, ymin=ymin.gzd, ymax=ymax.gzd), position=position_dodge(width=.2)) + geom_hline(aes(y=0, color='dummy'))
# ggplot(df.fancy.diffs, aes(x=lpos3c, color=comparison, group=comparison)) + geom_pointrange(aes(y=y.refix, ymin=ymin.refix, ymax=ymax.refix), position=position_dodge(width=.2)) + geom_hline(aes(y=0, color='dummy'))
# 
# ggplot(df.plot, aes(lpos3c, gzd3c, group=cond, color=cond)) + stat_summary(fun.data="mean_se")

# ## plotting
# palette <- c("#D95F02","#7570B3","#1B9E77")
# gzd.ylim <- c(200, 440)
# refix.ylim <- c(0, 1)
# ebwidth = .25

# ### yos plotting
# ggplot(subset(df.plot, expt=="yos"), aes(lpos3c, gzd3c.resid, group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(aes(group=cond, color=cond), fun.y=mean, geom="path") + scale_y_continuous("Gaze duration") + scale_x_continuous("Landing position") + coord_cartesian(ylim=gzd.ylim) + scale_color_manual(values=palette)
# 
# ggplot(subset(df.plot, expt=="yos"), aes(lpos3c, gzd3c, group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(aes(group=cond, color=cond), fun.y=mean, geom="path") + scale_y_continuous("Gaze duration") + scale_x_continuous("Landing position") + coord_cartesian(ylim=gzd.ylim) + scale_color_manual(values=palette)
# ggsave("gzd_yos.pdf", height=2.5,width=4)
# 
# ggplot(subset(df.plot, expt=="yos"), aes(lpos3c, as.numeric(refix3c), group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(fun.y=mean, geom="line") + scale_y_continuous("Refixation rate") + scale_x_continuous("Landing position") + coord_cartesian(ylim=refix.ylim) + scale_color_manual(values=palette)
# ggsave("refix_yos.pdf", height=2.5,width=4)
# 
# ### yosdos plotting
# ggplot(subset(df.plot, expt=="yd"), aes(lpos3c, gzd3c, group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(aes(group=cond, color=cond), fun.y=mean, geom="path") + scale_y_continuous("Gaze duration") + scale_x_continuous("Landing position") + coord_cartesian(ylim=gzd.ylim) + scale_color_manual(values=palette)
# ggsave("gzd_yd.pdf", height=2.5,width=4)
# 
# ggplot(subset(df.plot, expt=="yd"), aes(lpos3c, as.numeric(refix3c), group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(fun.y=mean, geom="line") + scale_y_continuous("Refixation rate") + scale_x_continuous("Landing position") + coord_cartesian(ylim=refix.ylim) + scale_color_manual(values=palette)
# ggsave("refix_yd.pdf", height=2.5,width=4)
# 
# ### yosdosreplication plotting
# ggplot(subset(df.plot, expt=="ydr"), aes(lpos3c, gzd3c, group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(aes(group=cond, color=cond), fun.y=mean, geom="path") + scale_y_continuous("Gaze duration") + scale_x_continuous("Landing position") + coord_cartesian(ylim=gzd.ylim) + scale_color_manual(values=palette)
# ggsave("gzd_ydr.pdf", height=2.5,width=4)
# 
# ggplot(subset(df.plot, expt=="ydr"), aes(lpos3c, as.numeric(refix3c), group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(fun.y=mean, geom="line") + scale_y_continuous("Refixation rate") + scale_x_continuous("Landing position") + coord_cartesian(ylim=refix.ylim) + scale_color_manual(values=palette)
# ggsave("refix_ydr.pdf", height=2.5,width=4)
# 
# ### yd + ydr plotting
# ggplot(subset(df.plot, expt %in% c("yd", "ydr")), aes(lpos3c, gzd3c, group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(aes(group=cond, color=cond), fun.y=mean, geom="path") + scale_y_continuous("Gaze duration") + scale_x_continuous("Landing position") + coord_cartesian(ylim=gzd.ylim) + scale_color_manual(values=palette)
# ggsave("gzd_yd_ydr.pdf", height=2.5,width=4)
# 
# ggplot(subset(df.plot, expt %in% c("yd", "ydr")), aes(lpos3c, as.numeric(refix3c), group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(fun.y=mean, geom="line") + scale_y_continuous("Refixation rate") + scale_x_continuous("Landing position") + coord_cartesian(ylim=refix.ylim) + scale_color_manual(values=palette)
# ggsave("refix_yd_ydr.pdf", height=2.5,width=4)
# 
# ### all three expts together plotting
# ggplot(subset(df.plot, !(expt=="yos" & (lpos3c > 4))), aes(lpos3c, gzd3c, group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(aes(group=cond, color=cond), fun.y=mean, geom="path") + scale_y_continuous("Gaze duration") + scale_x_continuous("Landing position") + coord_cartesian(ylim=gzd.ylim) + scale_color_manual(values=palette)
# ggsave("gzd_all.pdf", height=2.5,width=4)
# 
# ggplot(subset(df.plot, !(expt=="yos" & (lpos3c > 4))), aes(lpos3c, gzd3c.resid, group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(aes(group=cond, color=cond), fun.y=mean, geom="path") + scale_y_continuous("Gaze duration") + scale_x_continuous("Landing position") + coord_cartesian(ylim=gzd.ylim) + scale_color_manual(values=palette)
# ggsave("gzd_all_resid.pdf", height=2.5,width=4)
# 
# ggplot(subset(df.plot, !(expt=="yos" & (lpos3c > 4))), aes(lpos3c, as.numeric(refix3c), group=cond, color=cond)) + stat_summary(fun.data="mean_se", geom="errorbar", width=ebwidth) + stat_summary(fun.y=mean, geom="line") + scale_y_continuous("Refixation rate") + scale_x_continuous("Landing position") + coord_cartesian(ylim=refix.ylim) + scale_color_manual(values=palette)
# ggsave("refix_all.pdf", height=2.5,width=4)
