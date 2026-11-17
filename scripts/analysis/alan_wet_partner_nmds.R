# Each map uses one common source (one source for each edge for predator maps).
library(vegan)
root <- "results/all_species_forest/wet_partner_composition"
links <- read.csv(gzfile(file.path(root,"links.csv.gz")),stringsAsFactors=FALSE)
sources <- read.csv(file.path(root,"sources.csv"),stringsAsFactors=FALSE)
w <- read.csv(file.path(root,"species.csv"),stringsAsFactors=FALSE)
sources <- subset(sources,candidates>=2 & references>=2 & partners>=3 & plants>=6)
set.seed(9112)
coords <- list()
diagnostics <- list()
png(file.path(root,"nmds.png"),width=1950,height=750,res=150)
par(mfrow=c(1,3),mar=c(3.8,3.8,3.5,.8))
for(layer in unique(links$layer)){
    pool <- sources[sources$layer==layer,]
    pool <- pool[order(-pool$plants,-pool$partners,pool$flower_source,pool$diet_source),]
    plotted <- FALSE
    for(i in seq_len(nrow(pool))){
        s <- pool[i,]
        g <- subset(links,links$layer==s$layer & flower_source==s$flower_source & diet_source==s$diet_source)
        m <- table(g$sp,g$partner)>0
        if(nrow(unique(m))<4) next
        d <- vegdist(m,method="jaccard",binary=TRUE)
        if(length(unique(d))<3) next
        fit <- metaMDS(d,k=2,weakties=FALSE,try=5,trymax=15,trace=FALSE,wascores=FALSE,autotransform=FALSE)
        xy <- as.data.frame(scores(fit,display="sites"))
        xy$sp <- rownames(xy)
        xy$candidate <- w$candidate[match(xy$sp,w$sp)]=="True"
        xy$layer <- layer
        xy$flower_source <- s$flower_source
        xy$diet_source <- s$diet_source
        coords[[length(coords)+1]] <- xy
        diagnostics[[length(diagnostics)+1]] <- cbind(s,stress=fit$stress,plotted=!plotted)
        if(!plotted){
            plot(xy[,1:2],type="n",asp=1,main=switch(layer,visitors="Flower visitors",predators="Predators via visitors",adult_prey_predators="Predators via adult prey"))
            points(xy[,1:2],pch=ifelse(xy$candidate,17,1),col=ifelse(xy$candidate,"#AF503C","#55707C"),cex=1)
            legend("bottomleft",c("Decline candidate","Comparison"),pch=c(17,1),col=c("#AF503C","#55707C"),bty="n",cex=.7)
            mtext(sprintf("%d plants  |  stress %.2f",nrow(xy),fit$stress),side=3,line=.3,cex=.7)
            plotted <- TRUE
        }
    }
}
dev.off()
write.csv(do.call(rbind,coords),file.path(root,"ordination.csv"),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(root,"ordination_sources.csv"),row.names=FALSE)
