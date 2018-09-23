source("https://bioconductor.org/biocLite.R")
biocLite("rCGH")
library(rCGH)

data<-hg19
data$p<- data$centromerStart
data$q<- data$length - data$centromerEnd

library(tidyr)

data2<- gather(data, "arm", "length.arm", 6:7)

names(data2)

data3<- data2 %>% 
  dplyr::arrange(chrom) %>%
  unite("chromosome_arm", "chrom", "arm", sep="")

data3$chromosome<- parse_number(data3$chromosome_arm)

import<-fread(file.choose())
import<-totRes

View(import[1:100,])
head(import)

resList<-list()
for (i in 1:nrow(import)){
  chr<- which(import$Chromosome[i] == data$chrom)[1]
  sel<- dplyr::filter(data, chrom == chr)
  arm<- if(import$Start.Position[i] <= sel$centromerStart) {"p"
  } else if (import$Start.Position[i] >= sel$centromerEnd) {"q"}
  id<- paste0(chr, arm)
  resList[[i]]<-id
}

res<-unlist(resList)
final<-cbind(import,res)


final$length<- final$End.Position - final$Start.Position

percList<-list()
for(i in 1:nrow(final)){
  chrarm<-final$res[i]
  p<- final$length[i] / data3$length.arm[data3$chromosome_arm==chrarm] 
  percList[[i]]<- p
}

final2<- cbind(final, round(unlist(percList),4))




hist(final2$V2,breaks=100)


explore<-dplyr::filter(final2, length > 1000000)

explore<- dplyr::filter(explore, Seg.CN > 0.145 | Seg.CN < -0.177)


hist(explore$V2,breaks = 100)
