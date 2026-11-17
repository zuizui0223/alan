# Run from the alan repository root. Source archive DOI 10.5281/zenodo.4540407.
library(lme4)
root <- "results/all_species_forest"
out <- file.path(root,"wet_environment")
load("../gbif_annual/giavi.RData")
DAT$plab.log.z <- as.numeric(scale(log(DAT$plab)))
fit <- lmer(log1p(dVis) ~ plant*treatment + insect_group + plab.log.z + (1|pair), data=DAT)
lookup <- read.csv(gzfile(file.path(root,"decline_traits/moisture_evidence.csv.gz")))
estimates <- lapply(levels(factor(DAT$plant)),function(label) {
  g <- DAT[DAT$plant==label,]
  new <- g[rep(1,2),,drop=FALSE]
  new$treatment <- factor(c("C","T"),levels=levels(factor(DAT$treatment)))
  x <- model.matrix(delete.response(terms(nobars(formula(fit)))),new,
                    xlev=.getXlevels(terms(nobars(formula(fit))),model.frame(fit)))
  diff <- x[2,names(fixef(fit))]-x[1,names(fixef(fit))]
  effect <- sum(diff*fixef(fit))
  se <- sqrt(drop(t(diff) %*% as.matrix(vcov(fit)) %*% diff))
  name <- gsub("_"," ",label,fixed=TRUE)
  i <- if(grepl("^[A-Z][a-z]+ [a-z-]+$",name) && !endsWith(name," sp")) match(name,lookup$species) else NA_integer_
  data.frame(label=label,screened=!is.na(i),candidate=lookup$candidate[i],moisture=lookup$indicator[i],
             control_meadow_pairs=length(unique(g$site[g$treatment=="C"])),lit_meadow_pairs=length(unique(g$site[g$treatment=="T"])),
             shared_meadow_pairs=length(intersect(g$site[g$treatment=="C"],g$site[g$treatment=="T"])),
             estimate=effect,wald_low=effect-1.96*se,wald_high=effect+1.96*se)
})
write.csv(do.call(rbind,estimates),file.path(out,"giavi_diurnal.csv"),row.names=FALSE,na="")
# Equal weight per meadow pair, retaining sampling rounds where the plant occurs on both sides.
focal <- subset(DAT,plant=="Cirsium_oleraceum")
v <- aggregate(dVis~site+repeate+treatment,data=focal,FUN=sum)
matched <- merge(subset(v,treatment=="T"),subset(v,treatment=="C"),by=c("site","repeate"),suffixes=c("_lit","_dark"))
matched$difference <- log1p(matched$dVis_lit)-log1p(matched$dVis_dark)
by_pair <- aggregate(difference~site,matched,mean)
by_pair$matched_rounds <- as.integer(table(factor(matched$site,levels=by_pair$site)))
write.csv(by_pair,file.path(out,"cirsium_meadow_pairs.csv"),row.names=FALSE)
signs <- expand.grid(rep(list(c(-1,1)),nrow(by_pair)))
permuted <- as.matrix(signs) %*% by_pair$difference/nrow(by_pair)
print(c(delta=mean(by_pair$difference),p=mean(abs(permuted)>=abs(mean(by_pair$difference))-1e-10)))
