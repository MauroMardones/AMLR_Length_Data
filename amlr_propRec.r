# from May 18 2021
# modified for July 2021 "AMLR_Krill_LFD_data.csv" column names
#path.nm <- "/users/noaa/GYM_2021/amlr_lf_data/"
path.nm <- "/users/noaa/2022_propRec/"
rec.size <- 27 # 27 if recruits < 36mm, 31 if < 40mm
dat.all <- # net length frequencies including zero hauls
           read.csv(file=paste("AMLR_Krill_LFD_data.csv",sep=""),sep=",",
	   header=T,stringsAsFactors=F)
  dat.all$length[dat.all$length==0] <- NA # replace 840 zero lengths
  dat.all$amount[dat.all$amount==0] <- NA # replaced 683 zero amounts
dat.l <- subset(dat.all,dat.all$leg=="A" | dat.all$leg =="D") # A and D legs only
dat.l <- subset(dat.l,dat.l$amlr_area=="EI" | dat.l$amlr_area =="SA" |
                dat.l$amlr_area=="WA" | dat.l$amlr_area =="JI")
dat.l <-cbind(Year=substr(dat.l$AMLR_Cruise,5,8),dat.l)
dat.l <- dat.l[-which(is.na(dat.l$amount)),]
################
sqkm.areas <- c(43865,18151,24479,38524)
names(sqkm.areas) <- names.area <- c("EI","JI","SA","WA")
sqkm.all <- sum(sqkm.areas)
################
# abundance scaled up from individuals measured to captured
haul.totadults <- tapply(dat.l$Krill_total,list(dat.l$AMLR_Station),
                  unique)
haul.len.meas <- tapply(dat.l$amount,                 # Numbers measured in each length-class
                   list(dat.l$AMLR_Station,dat.l$length),  # by sample.
                   sum,na.rm=T)                       # (dim = c(2492,51)) 
haul.meas <- tapply(dat.l$amount,                     # Total numbers measured for length
                   list(dat.l$AMLR_Station),          # (dim = 2492)
                   sum,na.rm=T)                         
haul.len.meas[is.na(haul.len.meas)] <- 0
haul.len.meas <- haul.len.meas[,2:ncol(haul.len.meas)] # remove -1 length
haul.scalor <- array(dim=dim(haul.meas))              # Number by which to scale up measured
                 names(haul.scalor) <-                # individuals in each sample.(dim = 2492)
                 names(haul.meas)
  for(i.id in 1:length(haul.meas))
    ifelse(haul.totadults[i.id]>0,
           haul.scalor[i.id] <- haul.totadults[i.id]/haul.meas[i.id],
           haul.scalor[i.id] <- NA
           )         
################
  # Length-weight relationship
# length-wt from Farber-Lorda, 1994
# wt_len <- 0.00000503*(12:60)^3.283
# length-wt from Hewitt et al 2004
wt_len = 0.000002236*(12:60)^3.314; names(wt_len) <- 12:60

################
# LENGTH totals, scaled up from measured individuals per haul
haul.len.scaled <- 
                   array(dim=dim(haul.len.meas))      # Scaled-up numbers of individuals by length
  dimnames(haul.len.scaled) <-                        # for each sample.(dim = c(2492,51))
                   dimnames(haul.len.meas)
  for(i.id in 1:length(haul.meas))
    haul.len.scaled[i.id,] <- haul.len.meas[i.id,] * haul.scalor[i.id]
  haul.len.scaled[is.na(haul.len.scaled)] <- 0  

region.id <- tapply(as.character(dat.l$amlr_area),dat.l$AMLR_Station,paste)
leg.id <- tapply(as.character(dat.l$leg),dat.l$AMLR_Station,paste)
  for(i.id in 1:length(region.id)){
    region.id[i.id] <- region.id[[i.id]][1]
    leg.id[i.id] <- leg.id[[i.id]][1]
    }
region.id <- as.character(region.id) # change from array to character
leg.id <- as.character(leg.id)
yr.id <- substr(rownames(haul.len.meas),5,8)
haul.len.meas <- 
               cbind(region.id,yr.id,leg.id,as.data.frame(haul.len.meas),
               stringsAsFactors=T)                      # (dim = c(2492,52))

haul.len.meas$fr <-apply(as.matrix(haul.len.meas[,4:rec.size]),1,sum) /   # rec.size 27 is 35mm krill
                   apply(as.matrix(haul.len.meas[,4:ncol(haul.len.meas)]),1,sum)
haul.len.scaled <- 
               cbind(region.id,yr.id,leg.id,as.data.frame(haul.len.scaled),
               stringsAsFactors=T)                      # (dim = c(2492,49))
haul.len.scaled <- haul.len.scaled[,-c(3,4)] # strip out 2mm, 4mm length columns

haul.len.scaled$fr <-apply(as.matrix(haul.len.scaled[,4:rec.size]),1,sum) /   # rec.size 27 is 35mm krill
                   apply(as.matrix(haul.len.scaled[,4:ncol(haul.len.scaled)]),1,sum)


mean(haul.len.meas$fr,na.rm=TRUE)
# [1] 0.2152344 (when rec.size = 27)
sd(haul.len.meas$fr,na.rm=TRUE)
# [1] 0.3181458 

# Scaled mean by volume used (sensu D.K.)

mean(haul.len.scaled$fr,na.rm=TRUE)
# [1] 0.2522668 (when rec.size = 27)
sd(haul.len.scaled$fr,na.rm=TRUE)
# [1] 0.3384757

save("haul.len.meas", file = "meanlenghtKrillAMLR.RData")
save("haul.len.scaled", file = "meanscaledKrillAMLR.RData")




