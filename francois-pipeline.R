library(popcycle)
library(plotrix)

#convert -delay 0.1 -loop 0 -monitor $(ls *.class.gif | sort -n | perl -ne 'chomp($_); print $_ . " "') animation.gif

list <- list.files("~/Desktop/VCcal/", pattern='.evt', recursive=T, full.names=T)

png("~/Desktop/VCbeadsD2.png",width=20, height=20, unit='in', res=200)
par(mfrow=c(5,5), mar=c(2,1,1,1), oma=c(1,1,0,0))
for(file in list){
  #file <- list[2]
  print(file)
  evt <- readSeaflow(file, transform=T)
  plot.cytogram(evt, "fsc_small", "D2",  main=paste(basename(file)))
}
dev.off()

#####################
# 1. OPP FILTRATION #
#####################
cruise <- "MBARI_1"
path <- "/Volumes/data/data/seaflow/refilter/"
db <- paste0(path,cruise, "/",cruise,".db")
#make.popcycle.db(db)  # Create a popcycle SQLite3 database file

# opp.table <- get.opp.table(db)
# opp.table$date <- as.POSIXct(opp.table$date, format = "%FT%T", tz = "GMT")
# opp.table2 <- subset(opp.table, filter_id=="b6b02676-e6db-460e-a788-7aa26ebad1d8")
# plot(opp.table2$date, 100*opp.table2$evt_count/opp.table2$all_count)
# plot(opp.table2$date, 100*opp.table2$opp_evt_ratio)


#################
# 2. UPLOAD SFL #
#################
for(i in c(9,10,14,15,17)){
  #i <- 17
  cruise <- paste0("SCOPE_",i)
  print(cruise)
  path <- "/Volumes/data/data/seaflow/refilter/"
  path.sfl <- paste0("~/Desktop/scope/",cruise,"/cruise.sfl")
  db <- paste0(path,cruise, "/",cruise,".db")
  
  df <- read.delim(path.sfl)
  gga <- FALSE
  if(mean(df$LAT, na.rm=T) > 1000) gga <- TRUE
  save.sfl(db=db, cruise=cruise, sfl.file=path.sfl, gga = gga)
  sfl <- get.sfl.table(db)
  print(sfl[1,])
}



########################
# 3. SET GATING PARAMS #
########################
# scope_16 wrong Year day
# scope_16 regate subarctic
cruise <- "SCOPE_16"
path <- "/Volumes/data/data/seaflow/refilter-2016-10-11/"
opp.dir <- paste0(path,cruise, "/",cruise,"_opp")
vct.dir <- paste0(path,cruise, "/", cruise, "_vct")
db <- paste0(path,cruise, "/",cruise,".db")

opp.list <- get.opp.files(db)

freq <- round(seq(75, length(opp.list), length.out=24))

OPP <- NULL
par(mfrow=c(4,6))

for(i in freq){
  opp.name <- opp.list[i]
  try(opp <- get.opp.by.file(opp.dir, opp.name))#, vct.dir=vct.dir, pop= "prochloro"))
  try(plot.cytogram(opp, "fsc_small","chl_small", main=paste(basename(opp.name))))
  OPP <- rbind(OPP, opp)
}

gates.log <- add.manual.classification(OPP, "beads", "fsc_small", "pe")
gates.log <- add.manual.classification(OPP, "synecho", "fsc_small", "pe", gates.log)
#gates.log <- add.manual.classification(OPP, "croco", "fsc_small", "pe", gates.log)
#gates.log <- add.auto.classification("synecho", "fsc_small", "pe", position=c(FALSE,TRUE), gates=c(2.0,NA), scale=0.975, min.pe=7, gates.log=gates.log)
#gates.log <- add.auto.classification("prochloro", "fsc_small", "chl_small",  position=c(FALSE,TRUE), gates=c(2.0,0.5), scale=0.975, gates.log=gates.log)
gates.log <- add.manual.classification(OPP, "prochloro", "fsc_small", "chl_small", gates.log)
gates.log <- add.manual.classification(OPP, "picoeuk", "fsc_small", "chl_small", gates.log)



### Check classification on a subset
par(mfrow=c(4,6))
for(i in freq){
  opp.name <- opp.list[i]
  opp <- try(get.opp.by.file(opp.dir, opp.name))
  opp <- try(classify.opp(opp, gates.log))
  try(plot.vct.cytogram(opp, para.x="fsc_small", para.y="chl_small", main=paste(basename(opp.name))))
}

### CLASSIFY
gating.id <- save.gating.params(db, gates.log)
classify.opp.files(db, cruise, opp.dir, opp.list[c(1:length(opp.list))], vct.dir, gating.id=gating.id); gating.id <- NULL



#######################
### 4. Check gating ###
#######################
cruise <- "SCOPE_19"
path <- "/Volumes/data/data/seaflow/refilter-2016-10-11/"
opp.dir <- paste0(path,cruise, "/",cruise,"_opp")
vct.dir <- paste0(path,cruise, "/", cruise, "_vct")
db <- paste0(path,cruise, "/",cruise,".db")
opp.list <- get.opp.files(db)

freq <- round(seq(25, length(opp.list), length.out=24))
par(mfrow=c(4,6))
for(i in freq){
  opp.name <- opp.list[i]
  try(plot.vct.cytogram.by.file(opp.dir,  vct.dir, opp.name, para.y="chl_small", main=paste(basename(opp.name))))
  #try(plot.opp.cytogram.by.file(opp.dir,  opp.name, main=paste(opp.name), para.y="chl_small"))
}

stat <- get.stat.table(db, flag=F)
par(mfrow=c(1,1))
plot.time(stat, popname='prochloro', param='abundance')



### CHECK TABLES IN POPCYCLE.DB
# sqlite3 path/to/mydb.db
# drop table vct;
# sqlite3 path/to/mydb.db <path/to/popcycle.sql
for(i in c(9,10,14,15,17)){
  
  #i <- 17
  cruise <- paste0("SCOPE_",i)
  print(cruise)
  path <- "/Volumes/data/data/seaflow/refilter/"
  db <- paste0(path,cruise, "/",cruise,".db")
  
  sql <- paste0("SELECT * FROM vct")
  #sql <- paste0("drop table sfl")
  #sql <- paste0("drop table vct")
  #sql <- paste0("drop view stat")
  
  t <- sql.dbGetQuery(db, sql)
  print(t[1,2])
  system(paste("sqlite3", db, "<~/Documents/DATA/Codes/popcycle/inst/sql/popcycle.sql"))
  
}





###################
### 5. OUTLIERS ###
###################
library(popcycle)

cruise <- "SCOPE_16"
path <- "/Volumes/data/data/seaflow/refilter/"
db <- paste0(path,cruise, "/",cruise,".db")
stat <- get.stat.table(db, flag=F)
#stat <- subset(stat, flag == 0)
plot.time(stat, popname='prochloro', param='fsc_small_mean')
#  write.csv(stat, "~/Desktop/stat.csv", quote=F, row.names=F)
stat$time <- as.POSIXct(stat$time, format = "%FT%T", tz = "GMT")
stat$flag <- 0
sfl <- get.sfl.table(db)
sfl$date <- as.POSIXct(sfl$date, format = "%FT%T", tz = "GMT")


par(mfrow=c(1,1), pty="m")


### remove FLOW_RATE outliers
para <- "flow_rate"
plot(sfl[,"date"], sfl[,para],type="p", ylab="Flow Rate", xlab="time", main='FLOW RATE outliers ')
fact.sd <- 0.5
model <- smooth.spline(sfl[,"date"], sfl[,para], spar=1.5)
lines(sfl[,"date"], fitted(model),lty=2, col=3)
res <- residuals(model)
pre.out <- which(res < -fact.sd*sd(res) | res > fact.sd*sd(res))
id <- sfl[pre.out,"file"]; out <- which(!is.na(match(stat[,"file"], id)))
points(sfl[pre.out,"date"], sfl[pre.out,para], col=2)

stat[out,'flag'] <- 1

### remove EVENT_RATE outliers
coinc <- 18000
para <- "event_rate"
plot(sfl[,"date"], sfl[,para],type="p", ylab="Event Rate", xlab="time", main='EVENT_RATE outliers ')
abline(h=coinc, col=3, lty=2)
pre.out <- which(sfl[,para] > coinc)
points(sfl[pre.out,"date"], sfl[pre.out,para], col=2)
id <- sfl[pre.out,"file"]; out <- which(!is.na(match(stat[,"file"], id)))

stat[out,'flag'] <- 1

### remove OPP FILTRATION outliers
df <- subset(stat,flag==0)
para <- "opp_evt_ratio"
plot(df[,"time"], df[,para],type="p", ylab="OPP/EVT ratio", xlab="time", main='OPP FILTRATION outliers ')
fact.sd <- 2.5
model <- smooth.spline(df[,"time"], df[,para], spar=0.125)
lines(df[,"time"], fitted(model),lty=2, col=3)
res <- residuals(model)
out <- which(res < -fact.sd*sd(res) | res > fact.sd*sd(res))
#out <- which(df[,para] < 0.005)
points(df[out,"time"], df[out,para], col=2)

stat[rownames(df[out,]), 'flag'] <- 2


# ## remove BEADS outliers
# print('remove BEADS outliers ')
#
# beads <- subset(stat, flag ==0 & pop=='beads')
# para <- beads[,"fsc_small_mean"]
# plot(as.POSIXct(beads[,"time"]), para,type="p", ylab="BEADS fsc_small", xlab="time", log='y')
# fact.sd <- 2
# model <- smooth.spline(beads[,"time"], para, spar=0.75)
# lines(beads[,"time"], fitted(model),lty=2,col='green')
# res <- residuals(model)
# pre.out <- which(res < -fact.sd*sd(res) | res > fact.sd*sd(res))
# points(beads[pre.out,"time"], beads[pre.out,"fsc_small_mean"],col=2)
# id <- beads[pre.out,"file"]; out <- which(!is.na(match(stat[,"file"], id)))
#
# stat[out, 'flag'] <- 2


### remove GATING outliers
lim.t <- c(min(stat$time,na.rm=T), max(stat$time,na.rm=T))
phyto <- c('prochloro', 'synecho')
parameters <- c("abundance", "fsc_small_mean")

par(mfrow=c(length(phyto), length(parameters)))

id <- NULL
for(para in parameters){
  print(para)
  for(i in phyto){
    p <- subset(stat, flag==0 & pop == i)
    plot(p$time, p[,para], xlim=lim.t, ylab=paste(para), xlab="time", main=paste(i))
    
    fact.sd <- 3.5
    model <- smooth.spline(p[,"time"], p[,para], spar=0.19)
    lines(p[,"time"], fitted(model),lty=2,col='green')
    res <- residuals(model)
    pre.out <- which(res < -fact.sd*sd(res) | res > fact.sd*sd(res))
    #pre.out <- which(p[,para] < 8)
    points(p[pre.out,"time"], p[pre.out,para],col='red')
    out <-  as.vector(unlist(data.frame(filename=unique(p[pre.out,"file"]))))
    id <- c(id, out)
    id <- unique(id)
  }
}


out <- which(!is.na(match(stat[,"file"], id)))
stat[out, 'flag'] <- 3


### plot CLEANED stat file
clean <- subset(stat, flag==0)
lim.t <- c(min(as.POSIXct(stat$time),na.rm=T), max(as.POSIXct(stat$time),na.rm=T))
phyto <- unique(stat$pop)

cex <- 1
par(mfrow=c(ceiling(length(phyto)),2), cex=cex, mar=c(2,4,2,3), oma=c(2,1,1,1))
n <- 1
for(i in phyto){
  p <- subset(clean, pop == i)
  print(i)
  if(nrow(p) > 0){
    plot(p$time, p$abundance, xlim=lim.t,xlab=NA,ylab=NA, main=paste(i),las=1)
    mtext(substitute(paste("Abundances (10"^{6},"cells L"^{-1},")")), 2, line=3, cex=cex)
    plot(p$time, p$fsc_small_mean, xlim=lim.t,xlab=NA,ylab=NA, main=paste(i),las=1)
    mtext("Light scattering", 2, line=3, cex=cex)
    if(nrow(p) < 10) mtext(paste("only ", nrow(p), "data points"), side=1, line=-4)
    if(n == length(phyto)) mtext("Time", side=1, line=3,cex=cex)
    n <- n + 1
  }
}




#### SAVE STAT TABLE AND SAVE OUTLIERS TO DB
# write.csv(stat, paste0(path,cruise, "/stat.csv"), quote=FALSE, row.names=FALSE)

df <- stat[match(unique(stat$file),stat$file),]
outlier <- data.frame(file=df$file, flag=df$flag)
reset.outlier.table(db)
save.outliers(db, cruise, table.name=outlier)

stat <- subset(stat, flag ==0 & pop != 'beads')
write.csv(stat[,c(3:5,9,11:14)], "~/Desktop/stat.csv", quote=FALSE, row.names=FALSE)
