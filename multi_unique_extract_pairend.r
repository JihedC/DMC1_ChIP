# This R script is used to extract multi-unique reads.
# revised on 20/05/2016
# revised on 12/10/2017
# args[1] is the summary stats file for multi_fq10.sam lines
# args[2] is the input sam file multi_fq10.sam
# args[3] is halfway saved the .RData file for multi-unique
# args[4] is the final output sam file without header.
# args[5] skip number, is the total lines of your header in sam file
# args[6] number of lines each time you want to read in R(this is trying to avoid big data problem.)
# args[7] total columns in sam file
# args[8] columns for the tag, 
# note: for single end, args[7] samfile has 20 columns
# for pair-end, args[7] is 21. 
# args[8] for single end is 15; pair-end is 16 instead.
args<-commandArgs(TRUE)
dim1 <- ceiling(read.table(args[1], header=F)[1,1]/as.numeric(args[6]))
infile <- args[2]
dat_out <- NULL
for(s in 1:dim1){
  print(s)
  dat_in <-  matrix(scan(infile, what="character",skip=as.numeric(args[5])+(s-1)*as.numeric(args[6]), nlines=as.numeric(args[6])), ncol=as.numeric(args[7]), byrow=T)
  colnames(dat_in) <- c("qname", "flag", "chr", "locs", "mapq", paste("tag",c(1:as.numeric(args[8])), sep=""))
  qnames.len <- length(unique(dat_in[,1]))
  test <-  by(dat_in, list(qname=dat_in[,1]), function(x){
    y <- subset(x, select= -qname)
    y})
  testfn <- function(ldat, ind){
    newd <- ldat[[ind]][1,]
    find <- as.matrix(newd,ncol=(as.numeric(args[7])-1))
  }
  ntest <- sapply(c(1:qnames.len), function(x) testfn(test,x))
  finaltest <- t(ntest)
  newnames <- names(test)
  totaldat <- cbind(newnames, finaltest)
  dat_out <- rbind(dat_out, totaldat)
  dat_out
}
colnames(dat_out) <- c("qname", "flag", "chr", "locs", "mapq", paste("tag",c(1:as.numeric(args[8])), sep=""))
rpi <- dat_out
save(rpi, file=args[3])
check.len0 <- dim(rpi)[1]
check.len1 <- length(unique(rpi[,1]))
s0 <- table(rpi[,1])
multi <- names(which(s0>1))
dindex <- rpi[,1]%in%multi
subdat <- rpi[dindex==T,]
olddat <- rpi[dindex==F,]
qnames.len <- length(unique(subdat[,1]))
test <-  by(subdat, list(qname=subdat[,1]), function(x){
  y <- subset(x, select= -qname)
  y})
ntest <- sapply(c(1:qnames.len), function(x) testfn(test,x))
subtest <- t(ntest)
newnames <- names(test)
newsubdat <- cbind(newnames, subtest)
colnames(newsubdat) <- c("qname", "flag", "chr", "locs", "mapq", paste("tag",c(1:as.numeric(args[8])), sep=""))
final.rpi <- rbind(olddat,newsubdat)
write.table(final.rpi, file = args[4], row.names = F, col.names=F, quote=F, sep="\t")
