library(Hmisc)
library(dplyr)

source("et-admin/et-admin.R")
source("etAnalyze/R/etAnalyze.R")

subj_excl <- readRDS("subj_excl.rds")
df_exclusions <- readRDS("df_exclusions.rds")

expt.yos <- get.expt.from.da1list("absDA1s.lst", "../yosPre.cnt",
                                  "Yosemite/Analysis/DataAbs/",
                                  min.fix = 80,
                                  ignore.num.fixes=T)
expt.yd <- get.expt.from.da1list("absDA1s.lst", "preCount.cnt",
                                 "YosDos/Analysis/",
                                 min.fix = 80,
                                 ignore.num.fixes=T)
expt.ydr <- get.expt.from.da1list("absDA1s.lst", "../preCount.cnt",
                                  "YosDosReplication/Analysis/DataAbs/",
                                  min.fix = 80,
                                  ignore.num.fixes=T)

df.second.yos <- second.saccade.info(expt.yos, 2) %>%
  mutate(expt = 'yos')
df.second.yd <- second.saccade.info(expt.yd, 2) %>%
  mutate(expt = 'yd')
df.second.ydr <- second.saccade.info(expt.ydr, 2) %>%
  mutate(expt = 'ydr')

df.second <- bind_rows(df.second.yos, df.second.yd, df.second.ydr) %>%
  left_join(df_exclusions, by=c("subj", "item", "expt")) %>%
  filter(blinks.good == T) %>% # exclude blinks only (not slow display changes)
  left_join(subj_excl, by = c("subj", "expt")) %>%
  filter(to_kick == F) %>% # exclude subjects who are excluded in main analysis
  select(-c(to_kick)) %>%
  mutate(subj = paste(expt, subj), # make subjects unique
         distance = x.next - x)

collapse.all.equals <- function(x) {
  if (length(x) == 0) {
    return(x)
  }
  stopifnot(all(x == x[1]))
  return(x[1])
}

## figure out pretarget lengths
df.pretarget.length <- df.second %>%
  filter(region.next == 3) %>% # we can use lpos on the target to make inferences about pretarget length
  mutate(pretarget.length = lpos + distance - lpos.next - 1) %>%
  group_by(item) %>%
  summarize(pretarget.length = collapse.all.equals(pretarget.length)) %>%
  ungroup()

stopifnot(nrow(df.pretarget.length) == 160)

## add in pretarget lengths
df.second <- df.second %>%
  left_join(df.pretarget.length, by=c("item")) %>%
  mutate(lposrel = distance + lpos,
         lpostarget = lposrel - pretarget.length - 1,
         launch = lpos - pretarget.length - 1)

saveRDS(df.second, file = "df.second.rds")
saveRDS(df.pretarget.length, "df.pretarget.length.rds")
