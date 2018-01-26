library(Hmisc)
library(dplyr)

source("et-admin/et-admin.R")

timing_max <- 12

cat("Yos timing exclusions\n")
df_dispchange_yos <- get.timing.exclusions("timing.txt", timing_max, 'timingAllList.lst', "Yosemite/Analysis/") %>%
  mutate(expt = 'yos')
cat("YD timing exclusions\n")
df_dispchange_yd <- get.timing.exclusions("timing.txt", timing_max, 'timingAllList.lst', "YosDos/Analysis/") %>%
  mutate(expt = 'yd')
cat("YDR timing exclusions\n")
df_dispchange_ydr <- get.timing.exclusions("timing.txt", timing_max, 'timingAllList.lst', "YosDosReplication/Analysis/") %>%
  mutate(expt = 'ydr')
df_dispchange <- bind_rows(df_dispchange_yos, df_dispchange_yd, df_dispchange_ydr)

df_blinks_yos <- get.blink.exclusions("absDA1s.lst", 160, "Yosemite/Analysis/DataAbs/") %>%
  mutate(expt = 'yos')
df_blinks_yd <- get.blink.exclusions("absDA1s.lst", 160, "YosDos/Analysis/") %>%
  mutate(expt = 'yd')
df_blinks_ydr <- get.blink.exclusions("absDA1s.lst", 160, "YosDosReplication/Analysis/DataAbs/") %>%
  mutate(expt = 'ydr')
df_blinks <- bind_rows(df_blinks_yos, df_blinks_yd, df_blinks_ydr)

### ensure that blink and timing labels align
df_dispchange <- df_dispchange %>%
  mutate(label = gsub(".asc", "", label))
df_blinks <- df_blinks %>%
  mutate(label = gsub(".da1", "", label))
stopifnot(identical(df_dispchange[c("subj", "item", "label")],
                    df_blinks[c("subj", "item", "label")]))

df_exclusions <- df_blinks %>%
  mutate(blinks.good = good) %>%
  select(-good) %>%
  left_join(df_dispchange, by = c("expt", "subj", "item", "label")) %>%
  mutate(dcgood = good) %>%
  select(-good) %>%
  mutate(good = blinks.good & dcgood) %>%
  select(subj, label, item, expt, blinks.good, dcgood, good)

# figure out subjects to exclude
subj_excl <- df_exclusions %>%
  group_by(expt, label, subj) %>%
  summarize(blinks.good = mean(blinks.good),
            dcgood = mean(dcgood),
            good = mean(good)) %>%
  ungroup() %>%
  mutate(to_kick = good < .5 | blinks.good < 2/3 | grepl("-", label)) %>%
  select(subj, expt, to_kick)

exclusion_summary <- subj_excl %>%
  group_by(expt) %>%
  summarize(subjs_kept = sum(!to_kick))
stopifnot(all(exclusion_summary$subjs_kept == 40)) # 40 subjs in each expt

saveRDS(subj_excl, "subj_excl.rds")
saveRDS(df_exclusions, "df_exclusions.rds")
