# Original author data, not a prior meta-analysis extraction.
# Source readme: site identifies meadow pair; repeate sampling event.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
source <- args[1]
out <- args[2]
stopifnot(unname(tools::md5sum(source))=='f60b21d82c163cbb53520fda0bc61985')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
env <- new.env()
load(source, envir=env)
d <- env$DAT
stopifnot(nrow(d)==1884, length(unique(d$site))==6,
          length(unique(d$plant))==21, all(d$treatment %in% c('C','T')),
          all(is.finite(d$dVis)), all(d$dVis>=0))
events <- unique(d[,c('site','repeate','treatment')])
stopifnot(nrow(events)==74, all(table(events$treatment)==37))
keys <- lapply(c('C','T'),function(t) sort(paste(events$site[events$treatment==t],events$repeate[events$treatment==t])))
stopifnot(identical(keys[[1]],keys[[2]]))
# Every observed plant-event has all three guild counts: no missing-guild zero fill.
stopifnot(all(table(paste(d$site,d$repeate,d$treatment,d$plant))==3))
event <- aggregate(dVis~site+repeate+treatment+insect_group,d,sum)
all <- aggregate(dVis~site+repeate+treatment,d,sum)
all$insect_group <- 'all_three_orders'
event <- rbind(event,all[,names(event)])
unit <- aggregate(dVis~site+treatment+insect_group,event,mean)
write.csv(unit,file.path(out,'giavi2021_pair_arm_means.csv'),row.names=FALSE)
write.csv(event,file.path(out,'giavi2021_sampling_events.csv'),row.names=FALSE)
endpoints <- sort(unique(unit$insect_group))
lit <- sapply(endpoints,function(k) {x<-subset(unit,treatment=='T' & insect_group==k);x$dVis[order(x$site)]})
dark <- sapply(endpoints,function(k) {x<-subset(unit,treatment=='C' & insect_group==k);x$dVis[order(x$site)]})
stopifnot(all(dim(lit)==c(6,4)),all(dim(dark)==c(6,4)),all(colMeans(lit)>0),all(colMeans(dark)>0))
# Delta method for ln(mean lit/mean dark), including paired meadow covariance.
n <- nrow(lit);ml<-colMeans(lit);md<-colMeans(dark)
v <- cov(lit)/n/outer(ml,ml)+cov(dark)/n/outer(md,md)-
     cov(lit,dark)/n/outer(ml,md)-cov(dark,lit)/n/outer(md,ml)
stopifnot(max(abs(v-t(v)))<1e-12,min(eigen(v,symmetric=TRUE)$values)> -1e-12)
effect <- data.frame(effect_id=paste0('Giavi2021_',endpoints),study_id='Giavi2021',
  programme_id='KnopGiavi_overlap_unresolved',experiment_id='Giavi2016_meadows',
  control_group_id='Giavi2016_dark6',endpoint=paste0('daytime_selected21_visits_',endpoints),
  independent_unit='meadow',paired_units=6,design='paired_meadows',
  mean_lit=ml,sd_lit=apply(lit,2,sd),n_lit=6,mean_dark=md,sd_dark=apply(dark,2,sd),n_dark=6,
  paired_arm_covariance=diag(cov(lit,dark)),lnRR=log(ml/md),sampling_variance=diag(v),
  measurement='mean visits per sampling event to author-selected 21 plant labels; 6 rounds/event',
  adjustment='raw visitation not adjusted for flower abundance',
  source='10.5281/zenodo.4540407 original data.Rdata DAT',
  verification='paired unit reconstruction; independent extraction pending',row.names=NULL)
write.csv(effect,file.path(out,'giavi2021_effects.csv'),row.names=FALSE)
grid<-expand.grid(i=seq_along(endpoints),j=seq_along(endpoints))
covariance<-data.frame(effect_id_1=effect$effect_id[grid$i],effect_id_2=effect$effect_id[grid$j],
  sampling_covariance=v[as.matrix(grid)],method='paired_multivariate_delta_same_six_meadow_pairs')
write.csv(covariance,file.path(out,'giavi2021_covariance.csv'),row.names=FALSE)
write.csv(data.frame(source_doi='10.5281/zenodo.4540407',file='data.Rdata',
  md5=unname(tools::md5sum(source)),raw_rows=1884,plant_labels=21,
  meadow_pairs=6,physical_meadows=12,paired_sampling_events=37),
  file.path(out,'giavi2021_source_manifest.csv'),row.names=FALSE)
print(effect[,c('endpoint','mean_lit','mean_dark','lnRR','sampling_variance')])
