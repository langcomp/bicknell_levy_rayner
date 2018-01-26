all : plotting.Rout analyze_mogs.Rout get_stats.Rout get_CIs.Rout

# all : plotting.Rout mogs.Rout stats.Rout
.PHONY: all clean

exclusions.Rout : exclusions.R
	Rscript $< > $@ 2> exclusions.Rerr
get_saccade_data.Rout : get_saccade_data.R exclusions.Rout
	Rscript $< > $@ 2> get_saccade_data.Rerr
get_data.Rout : get_data.R exclusions.Rout get_saccade_data.Rout
	Rscript $< > $@ 2> get_data.Rerr
get_mogs.Rout : get_mogs.R get_saccade_data.Rout
	Rscript $< > $@ 2> get_mogs.Rerr
analyze_mogs.Rout : analyze_mogs.R get_mogs.Rout
	Rscript $< > $@ 2> analyze_mogs.Rerr
get_stats.Rout : get_stats.R get_data.Rout
	Rscript $< > $@ 2> get_stats.Rerr
get_CIs.Rout : get_CIs.R get_stats.Rout
	Rscript $< > $@ 2> get_CIs.Rerr
plotting.Rout : plotting.R get_data.Rout
	Rscript $< > $@ 2> plotting.Rerr
clean :
	rm -f exclusions.Rout exclusions.Rerr subj_excl.rds df_exclusions.rds
	rm -f get_saccade_data.Rout get_saccade_data.Rerr df.second.rds df.pretarget.length.rds
	rm -f get_data.Rout get_data.Rerr df.lmer.rds
	rm -f get_mogs.Rout get_mogs.Rerr mogs.rds
	rm -f analyze_mogs.Rout analyze_mogs.Rerr
	rm -f get_stats.Rout get_stats.Rerr stats.RData
	rm -f get_CIs.Rout get_CIs.Rerr CIs.RData
	rm -f plotting.Rout plotting.Rerr gzd_nolaunch_full.pdf refix_nolaunch_full.pdf gzd_launch_full.pdf refix_launch_full.pdf
