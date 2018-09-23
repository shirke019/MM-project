source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

library(GenomicRanges)

data<-fread(file.choose())

#explore focal region cutoffs, default= 3Mb ( https://biology.stackexchange.com/questions/3407/what-is-a-focal-copy-number-variation )
data2<-data
data2$len<-data$End.Position - data$Start.Position
data2<-dplyr::filter(data2, len<=3000000)
data2<-dplyr::filter(data2, Seg.CN > 0.145 | Seg.CN < -0.177 )


hist(data2$len, xlim=c(0,1000000), breaks=100000)
library(ggplot2)
ggplot(data2, aes(len))+ xlim(c(0,1000000)) + geom_histogram(bins=1000) + facet_wrap(~Chromosome)
plot(table(data2$Chromosome))




prepdata<-data[,-c(5)]

prepdata$len<- prepdata$End.Position - prepdata$Start.Position
colnames(prepdata) <- c("id",'chr','start','end','score',"len")

#Call AMPS and DELS on log2ratio thresholds 
prepdata$score[prepdata$score< 0.145 & prepdata$score> -0.177]<- 0
prepdata$score[prepdata$score>=0.145]<- 1
prepdata$score[prepdata$score<= -0.177]<- -1

#filter out normal regions
prepdata<-dplyr::filter(prepdata, score!=0)

#filter out focal regions
#prepdata<-dplyr::filter(prepdata, len >= 3000000)

#LOOP PREPATATION for each patient
patientsList<-unique(prepdata$id) 

setdiff(unique(data2$Sample), patientsList)# NO ALTERATIONS PATIENTS

#LOOP
ampList<-list()
delList<-list()
for( i in patientsList){
  #filter for specific patient
  patient<-dplyr::filter(prepdata, id == i)
  #split by AMPS and DELS
  patientAMPS<-dplyr::filter(patient, score==1)
  patientDELS<-dplyr::filter(patient, score== -1)
  
  #____ DELS _____
  if(nrow(patientDELS)>0) {
  BED.DELS<-GRanges(patientDELS)
  BED.DELS
  #View(BED.DELS)
  BED2.DELS <- reduce(BED.DELS, min.gapwidth=3000000)

  IDX <- findOverlaps(BED.DELS, BED2.DELS)
  IDX2 <- IDX[which(!duplicated(subjectHits(IDX))),] #Just assign things once
  mcols(BED2.DELS)$id[subjectHits(IDX2)] <- mcols(BED.DELS)$id[queryHits(IDX2)]
  #View(BED2.DELS)
  delList[[i]]<- as.data.frame(BED2.DELS)
  }
  
  #____ AMPS _____
  if(nrow(patientAMPS)>0) {
  BED.AMPS<-GRanges(patientAMPS)
  BED.AMPS
  #View(BED.AMPS)
  BED2.AMPS <- reduce(BED.AMPS, min.gapwidth=3000000)

  IDX <- findOverlaps(BED.AMPS, BED2.AMPS)
  IDX2 <- IDX[which(!duplicated(subjectHits(IDX))),] #Just assign things once
  mcols(BED2.AMPS)$id[subjectHits(IDX2)] <- mcols(BED.AMPS)$id[queryHits(IDX2)]
  #View(BED2.AMPS)
  ampList[[i]]<- as.data.frame(BED2.AMPS)
  }
}

#RESULTS DATA FRAMES with CALLS in BED FORMAT
ampRes<-Reduce(rbind, ampList) 

delRes<-Reduce(rbind, delList)

ampRes$type<- "amplifcation"
delRes$type<- "deletion"

totRes<-rbind(ampRes, delRes)

#________________PART 2____________________

library(rCGH)

data<-hg19
data$p<- data$centromerStart
data$q<- data$length - data$centromerEnd


import<-totRes

head(import)

resList<-list()
for (i in 1:nrow(import)){
  chr<- which(import$seqnames[i] == data$chrom)[1]
  sel<- dplyr::filter(data, chrom == chr)
  arm<- if(import$start[i] <= sel$centromerStart) {"p"
  } else if (import$start[i] >= sel$centromerEnd) {"q"}
  id<- paste0(chr, arm)
  resList[[i]]<-id
}

res<-unlist(resList)
final<-cbind(import,res)


percList<-list()
for(i in 1:nrow(final)){
  chrarm<-final$res[i]
  p<- final$width[i] / data3$length.arm[data3$chromosome_arm==chrarm] 
  percList[[i]]<- p
}

final2<- cbind(final, perc=round(unlist(percList),4))




hist(final2$perc,breaks=100)


explore<-dplyr::filter(final2, width > 50000)

hist(explore$perc,breaks = 1000)

plot(density(explore$perc, bw=0.005))

