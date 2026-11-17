# Run from the alan repository. Distances compare positive habitat/trait annotations.
library(vegan)
root <- "results/all_species_forest"
out <- file.path(root,"seven_moth")
w <- read.csv(file.path(root,"wet_reproduction/species.csv"),stringsAsFactors=FALSE)
w$case <- w$candidate=="True"
w$seven <- w$case & w$documented_transport=="True"
h <- read.csv(gzfile(file.path(out,"habitat_optima.csv.gz")),stringsAsFactors=FALSE)
orig <- read.csv(gzfile(file.path(root,"decline_traits/moisture_evidence.csv.gz")))
w$region <- orig$site_region[match(w$sp,orig$sp)]
features <- c("Pollination.syndrome","Generative.reproduction.type","Reproduction.type")
supported <- w$sp %in% h$sp & apply(w[,features],1,function(x) all(!is.na(x) & nzchar(x)))
coverage <- w[,c("sp","species","case","seven")]
coverage$used <- supported
write.csv(coverage,file.path(out,"ordination_support.csv"),row.names=FALSE)
f <- w[supported,]
jaccard <- function(a,b) 1-length(intersect(a,b))/length(union(a,b))
hab <- split(h$optimum,h$sp)
tokens <- lapply(features,function(k) strsplit(f[[k]],", ",fixed=TRUE))
n <- nrow(f)
D <- matrix(0,n,n,dimnames=list(f$sp,f$sp))
H <- D
for(i in seq_len(n-1)) for(j in (i+1):n){
    H[i,j] <- H[j,i] <- jaccard(hab[[f$sp[i]]],hab[[f$sp[j]]])
    v <- sapply(tokens,function(lst) jaccard(lst[[i]],lst[[j]]))
    D[i,j] <- D[j,i] <- mean(c(H[i,j],v))
}
set.seed(714)
ord <- metaMDS(as.dist(D),k=2,try=12,trymax=55,trace=0,autotransform=FALSE,wascores=FALSE)
xy <- as.data.frame(scores(ord,display="sites"))
xy <- cbind(f[,c("sp","species","case","seven","region")],xy)
write.csv(xy,file.path(out,"ordination.csv"),row.names=FALSE)
pairs <- t(combn(which(f$seven),2))
write.csv(data.frame(first=f$species[pairs[,1]],second=f$species[pairs[,2]],habitat_distance=H[pairs],trait_habitat_distance=D[pairs]),file.path(out,"focal_distances.csv"),row.names=FALSE)
# Within the negative, insect-vector cohort compare the seven to sets of the same size
# with the same distribution of principal survey regions. The test describes this cohort.
test <- function(mat){
  pool <- which(f$case & grepl("insect-pollination",f$Pollination.syndrome))
  counts <- table(f$region[f$seven])
  meanD <- function(ix) mean(as.dist(mat[ix,ix]))
  observed <- meanD(which(f$seven))
  null <- replicate(9999,{
    idx <- unlist(lapply(names(counts),function(k){
      opts <- pool[f$region[pool]==as.numeric(k)]
      opts[sample.int(length(opts),counts[[k]])]
    }))
    meanD(idx)
  })
  data.frame(focal=observed,region_matched_mean=mean(null),lower_tail_p=(1+sum(null<=observed+1e-12))/10000,
             pool=length(pool),focal_n=sum(f$seven),regions=length(counts),stress=ord$stress)
}
tt <- rbind(habitat=test(H),combined=test(D))
write.csv(tt,file.path(out,"similarity_tests.csv"))
png(file.path(out,"habitat_nmds.png"),width=1800,height=1250,res=180)
par(mar=c(4,4,3,1))
plot(xy$NMDS1,xy$NMDS2,type="n",asp=1,xlab="NMDS 1",ylab="NMDS 2",main="Habitat optima, pollen vectors and reproduction")
points(xy$NMDS1,xy$NMDS2,pch=16,col=ifelse(xy$case,"#38877F","#C9D1D4"),cex=.8)
a <- xy[xy$seven,]
points(a$NMDS1,a$NMDS2,pch=21,bg="#EFAC5D",cex=1.6)
text(a$NMDS1,a$NMDS2,labels=seq_len(nrow(a)),cex=.75)
legend("bottomright",legend=paste(seq_len(nrow(a)),a$species),bty="o",bg="white",cex=.65)
legend("topleft",legend=c("Other screened wet species","Decline candidates","Seven pollen carriers' plants"),pch=16,col=c("#C9D1D4","#38877F","#EFAC5D"),bty="n",cex=.7)
dev.off()
print(tt)

# Do decline candidates as a whole occupy a tighter part of wet habitat space?
# Keep both the number of cases and the sampling stratum of each label fixed.
strata <- orig[match(f$sp,orig$sp),c("site_region","island_band","site_band")]
groups <- split(seq_len(n),interaction(strata,drop=TRUE))
both <- Filter(function(ix) length(unique(f$case[ix]))==2,groups)
observed <- unlist(lapply(both,function(ix) ix[f$case[ix]]),use.names=FALSE)
within_mean <- function(m,idx) mean(as.dist(m[idx,idx]))
set.seed(714)
null <- replicate(9999,{
  ix <- unlist(lapply(both,function(g) sample(g,sum(f$case[g]))),use.names=FALSE)
  c(within_mean(H,ix),within_mean(D,ix))
})
actual <- c(within_mean(H,observed),within_mean(D,observed))
p <- (1+rowSums(null<=actual+1e-12))/10000
write.csv(data.frame(distance=c("habitat","combined"),candidate_n=length(observed),
  reference_n=sum(lengths(both))-length(observed),strata=length(both),observed=actual,
  matched_null=rowMeans(null),p_lower=p,q=p.adjust(p,"BH")),file.path(out,"candidate_similarity.csv"),row.names=FALSE)
