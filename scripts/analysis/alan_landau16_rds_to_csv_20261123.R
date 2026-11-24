# Source-only conversion of pinned author original RDS; no model code packages.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('RDS source and output CSV paths required')
x <- readRDS(args[[1]])
if (!is.data.frame(x)) stop('Original source object is not data.frame')
if (!all(c('treat','week','flume','abund_rate') %in% names(x))) stop('Original data schema changed')
write.csv(x, file=args[[2]], row.names=FALSE, na='')
cat('SOURCE_RDS_READ rows',nrow(x),'cols',paste(names(x),collapse=','),'\n')
