## combined analysis of Yosemite, YosDos, and YosDosReplication
library(Hmisc)
library(dplyr)

source("et-admin/et-admin.R")
source("etAnalyze/R/etAnalyze.R")

subj_excl <- readRDS("subj_excl.rds")
df_exclusions <- readRDS("df_exclusions.rds")

df.pretarget.length <- readRDS("df.pretarget.length.rds")

get.all.data <- function(path,
                         has.sub, # has subgaze info?
                         expt) {
  ## get all the data
  df.ffix <- read.subj.by.item('FFIX.TXT', path)
  df.single <- read.subj.by.item('SINGLE.TXT', path)
  df.gzd <- read.subj.by.item('GZD.TXT', path)
  df.ro <- read.subj.by.item('REGOUT.TXT', path)
  df.nfix <- read.subj.by.item('NFIX.TXT', path)
  df.lpos <- read.subj.by.item('LPOS.TXT', path)
  df.launch <- read.subj.by.item('LAUNCH.TXT', path)
  
  df.ffix.pre <- read.subj.by.item('FFIXPRE.TXT', path)
  df.single.pre <- read.subj.by.item('SINGPRE.TXT', path)
  df.gzd.pre <- read.subj.by.item('GZDPRE.TXT', path)
  df.ro.pre <- read.subj.by.item('REGOPRE.TXT', path)
  df.nfix.pre <- read.subj.by.item('NFIXPRE.TXT', path)
  df.lpos.pre <- read.subj.by.item('LPOSPRE.TXT', path)
  df.launch.pre <- read.subj.by.item('LNCHPRE.TXT', path)
  
  if (has.sub) {
    df.gzd.sb <- read.subj.by.item('GZDSB.TXT', path)
    df.nfix.sb <- read.subj.by.item('NFIXSB.TXT', path)
  } else {
    df.gzd.sb <- df.gzd
    df.gzd.sb[, 4:ncol(df.gzd.sb)] <- 0
    df.nfix.sb <- df.nfix
    df.nfix.sb[, 4:ncol(df.nfix.sb)] <- 0
  }
  
  ## make sure they align
  check.alignment <- function(df) {
    identical(df.ffix[, c("subj", "item", "cond")],
              df[, c("subj", "item", "cond")])
  }
  stopifnot(check.alignment(df.single)
            & check.alignment(df.gzd)
            & check.alignment(df.ro)
            & check.alignment(df.nfix)
            & check.alignment(df.lpos)
            & check.alignment(df.launch)
            & check.alignment(df.ffix.pre)
            & check.alignment(df.single.pre)
            & check.alignment(df.gzd.pre)
            & check.alignment(df.ro.pre)
            & check.alignment(df.nfix.pre)
            & check.alignment(df.lpos.pre)
            & check.alignment(df.launch.pre)
            & check.alignment(df.gzd.sb)
            & check.alignment(df.nfix.sb))
  
  ## combine
  df.all <- df.ffix[, c("subj", "item", "cond")] %>%
    mutate(expt = expt)
  
  for (pre.label in c("", ".pre")) {
    for (measure.name in c("ffix", "single", "gzd", 'ro',
                           'nfix', 'lpos', 'launch')) {
      for (region.num in seq(2,4)) {
        source.df <- eval(parse(text=paste0("df.", measure.name, pre.label)))
        df.all[, paste0(measure.name, region.num, pre.label)] <-
          source.df[, region.num+3]
      }
    }
  }
  
  df.all$gzd3.sb <- df.gzd.sb$region2
  df.all$nfix3.sb <- df.nfix.sb$region2

  return(df.all)
}

df.yos <- get.all.data('Yosemite/Analysis/', F, 'yos')
df.yd <- get.all.data('YosDos/Analysis/', T, 'yd')
df.ydr <- get.all.data('YosDosReplication/Analysis/', T,
                       'ydr')

df.lmer <- bind_rows(df.yos, df.yd, df.ydr)

## add condition names
conds <- c("None", "Right", "None", "Left") # correspondence to numbers 1-4
df.lmer <- df.lmer %>%
  filter(cond != 0) %>%
  mutate(cond = factor(conds[cond],
                       levels = c("None", "Right", "Left")))

## create new variables
df.lmer <- df.lmer %>%
  mutate(gzd3c = as.integer(gzd3 - gzd3.sb),
         nfix3c = nfix3 - nfix3.sb,
         single3c = ifelse(nfix3c == 1, gzd3c, NA),
         ffix3c = ifelse(nfix3c == 0, NA, ffix3.pre), # only true for the data we analyze
         skip3c = (nfix3c == 0),
         skip2.pre = (nfix2.pre == 0),
         skip3.pre = (nfix3.pre == 0),
         refix3 = (nfix3 > 1),
         refix3.pre = (nfix3.pre > 1),
         refix3c = (nfix3c > 1),
         lpos3c = ifelse(cond == "None", lpos3.pre, NA),
         lpos3c = ifelse(cond == "Right", lpos3.pre - 3, lpos3c),
         lpos3c = ifelse(cond == "Left", lpos3.pre + 3, lpos3c),
         single2.pre = ifelse(skip2.pre | single2.pre == 0, NA, single2.pre),
         ffix2.pre = ifelse(skip2.pre | ffix2.pre == 0, NA, ffix2.pre),
         gzd2.pre = ifelse(skip2.pre | gzd2.pre == 0, NA, gzd2.pre))

stopifnot(with(df.lmer %>% filter(cond != "Left"),
               identical(gzd3c, gzd3) &
                 identical(nfix3c, nfix3)))

## When eyedry excludes trials for having too long of a duration,
## it just makes them a zero. This causes a problem for
## calculating gzd3c, so we exclude these cases across the board.
df.lmer <- df.lmer %>%
  filter(!(gzd3==0 & nfix3 > 0) &
           !(gzd3.sb==0 & nfix3.sb > 0))

## remove extraneous dependent measures that aren't meaningful anyway
df.lmer <- df.lmer %>%
  select(subj, item, cond, expt,
         ffix2.pre, single2.pre, gzd2.pre, ro2.pre, nfix2.pre, lpos2.pre, launch2.pre, launch3.pre,
         gzd3c, nfix3c, single3c, ffix3c, skip3c, skip2.pre, skip3.pre, refix3c, lpos3.pre, lpos3c,
         nfix3.sb, lpos3) # keep these two just for a sanity check later

df.lmer <- df.lmer %>%
  ## exclude subjects
  left_join(subj_excl, by=c("expt", "subj")) %>%
  filter(to_kick == F) %>%
  select(-c(to_kick)) %>%
  
  ## we're only analyzing cases where target is fixated both pre- and post-shift
  filter(!skip3c) %>%
  filter(skip3.pre == F) %>% # would have skipped target
  filter(launch3.pre != -1) %>% # b/c we exclude skips, these are problems (e.g., short fixations)
  filter(lpos3.pre != -1) %>% # would-be skips according to lpos
  mutate(throwoff = (((cond == "Right") & (lpos3.pre < 4)) | ((cond == "Left") & (lpos3.pre > 4)))) %>%
  filter(throwoff == F) %>%
  filter(lpos3.pre != 0) %>% # display change triggered without landing on word
  
  ## get rid of problematic cases
  left_join(df_exclusions, by=c("expt", "subj", "item"))

nrow(df.lmer %>% filter(!good)) / nrow(df.lmer)

df.lmer <- df.lmer %>%
  filter(good == T) # exclude blinks, long fixations, and display change problems

## sanity checks
stopifnot(with(df.lmer,
               (min(gzd2.pre, na.rm = TRUE) >= 80) &
                 (min(single3c, na.rm = TRUE) >= 80) &
                 (min(ffix2.pre, na.rm = TRUE) >= 80) &
                 (min(single2.pre, na.rm = TRUE) >= 80) &
                 (min(gzd3c) >= 80)
               ))
stopifnot(with(df.lmer,
               (cond=="None" & lpos3c == lpos3) |
                 (cond=="Right" & lpos3c == lpos3) |
                 (cond=="Left" & nfix3.sb == 0 & lpos3c == lpos3) |
                 (cond=="Left" & nfix3.sb > 0 & lpos3c > lpos3)))

## outlier removal
remove.outliers <- function(x, width=3) {
  x[x == 0] <- NA
  z <- scale(x)
  x[abs(z)>width] <- NA
  return(x)
}

df.lmer <- df.lmer %>%
  mutate(subj = paste(expt, subj)) %>%
  group_by(subj) %>%
  mutate(gzd3c = remove.outliers(gzd3c),
         single3c = remove.outliers(single3c),
         ffix2.pre = remove.outliers(ffix2.pre),
         ffix3c = remove.outliers(ffix3c),
         single2.pre = remove.outliers(single2.pre),
         gzd2.pre = remove.outliers(gzd2.pre)) %>%
  ungroup()

# how many gzd3c datapoints were excluded?
sum(is.na(df.lmer$gzd3c)) / nrow(df.lmer)

# add pretarget length
df.lmer <- df.lmer %>%
  left_join(df.pretarget.length, by=c("item")) %>%
  mutate(front.half.pretarget = ifelse(pretarget.length == 3, lpos2.pre > 1, NA),
         front.half.pretarget = ifelse(pretarget.length == 4, lpos2.pre > 2, front.half.pretarget))

# exclude very far launch sites with little data (this only affects the full dataset analysis. the subset already restricts launch site tightly)
old_size <- nrow(df.lmer)

table(df.lmer$launch3.pre) # exclude launch sites with fewer than 20 instances -> keeping 1-11
df.lmer <- df.lmer %>%
  filter(launch3.pre > 0, # this is a single case
         launch3.pre < 12)
1-nrow(df.lmer)/old_size # excluding 0.8%

# (optionally) analyze tight dataset to rule out mislocated refixations, etc.
# change full_dataset to TRUE this out to generate analyses / figures for the full dataset without these restrictions
full_dataset <- TRUE
if (!full_dataset) {
  cat("Dataset size before excluding multiple fixations of pretarget and fixations on back half of pretarget: ")
  cat(nrow(df.lmer), "\n")
  df.lmer <- df.lmer %>%
    filter(!skip2.pre,
           launch2.pre != -1, # b/c we exclude skips, these are problems (e.g., short fixations)
           ro2.pre == 0,
           nfix2.pre == 1,
           front.half.pretarget)

  cat("Dataset size after those exclusions: ")
  cat(nrow(df.lmer), '\n')
}

saveRDS(df.lmer, file = "df.lmer.rds")
