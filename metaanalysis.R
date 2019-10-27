library(dplyr)
library(tidyr)
library(ggplot2)
library(metafor)
library(meta)
library(dmetar)
library(grid)

classi <- read.csv("/home/chit/Desktop/Hiwi_BMC/Overlapstudies/classified.csv", stringsAsFactors = F)
overlap <- read.csv("/home/chit/Desktop/Hiwi_BMC/Overlapstudies/overlap.csv", stringsAsFactors = F)

#adjust the scale of the data to patients
overlap <- overlap %>%
  mutate(total=ifelse(scale=="eyes", total/2, total),
         case=ifelse(scale=="eyes", case/2, case)) %>%
  mutate(case=ifelse(scale=="percentage", total*(case/100), case)) %>%
  mutate(studylab=paste0(study,"_",autoimmune))

##meta-analysis
seTE <- c()
TE <- c()
lower <- c()
upper <- c()
validauto <- c()
classification <- c()
auto <- unique(overlap$autoimmune)
##for every autoimmune disease
for (x in auto){
filtered <- overlap %>%
  filter(autoimmune==x & ocular=="yes")

##remove only one case of manifestation
om <- unique(filtered$manifestation)
for (g in om){
  if (sum(filtered$manifestation==g)==1) {
    filtered <- filtered %>%
      filter(manifestation!=g)
  }
}

#if there is no usable data, skip the current iteration
if (dim(filtered)[1]==0){
  next
}

#meta-analysis for each autoimmune disease
meta <- metaprop(event = as.integer(case),
                 n= as.integer(total),
                 data = filtered,
                 studlab = studylab,
                 method="Inverse",
                 byvar = filtered$manifestation,
                 comb.random = T,
                 comb.fixed = F,
                 title=x)

png(sprintf("/home/chit/Desktop/Hiwi_BMC/Overlapstudies/%s.png", x), width = 1000, height = 1200)
forest(meta)
#grid.text(x, .5, 1-1/sizefactor, gp=gpar(cex=2))
dev.off()

seTE <- c(seTE, meta$seTE.random)
TE <- c(TE, meta$TE.random)
lower <- c(lower, meta$lower.random)
upper <- c(upper, meta$upper.random)
validauto <- c(validauto, x)

y <- classi[classi$ad==x,]$Classification
classification <- c(classification, y)
}

##meta-analysis for catergories
autoall <- data.frame("seTE"=meta:::backtransf(seTE, sm="PLOGIT"), 
                      "TE"=meta:::backtransf(TE, sm="PLOGIT"), 
                      "lower"=meta:::backtransf(lower, sm="PLOGIT"),
                      "upper"=meta:::backtransf(upper, sm="PLOGIT"),
                      validauto, classification)



meta_autoall <- metagen(TE=TE,
                       seTE=seTE,
                       lower = lower,
                       upper = upper,
                       data=autoall,
                       studlab = validauto,
                       byvar = autoall$classification,
                       sm="Prevalence")

png("/home/chit/Desktop/Hiwi_BMC/Overlapstudies/overall.png", width = 800, height = 700)
forest(meta_autoall,
       comb.random = T,
       comb.fixed = F,
       leftcols = "studlab",
       colgap.forest.left = "40mm",
       col.diamond.random = "lightblue",
       col.square = "darkblue",
       col.square.lines = "darkblue")
dev.off() ##saved


##statistic for non-ocular
nonocular <- overlap %>%
  filter(ocular=="no")

nonocular_auto <- unique(nonocular$autoimmune)
for (x in nonocular_auto){
  filtered <- nonocular %>%
    filter(autoimmune==x)
  
  ##remove only one case of manifestation
  om <- unique(filtered$organ)
  for (g in om){
    if (sum(filtered$organ==g)==1) {
      filtered <- filtered %>%
        filter(organ!=g)
    }
  }
  
  if (dim(filtered)[1]==0){
    next
  }
  
  #meta-analysis
  meta <- metaprop(event = as.integer(case),
                   n= as.integer(total),
                   data = filtered,
                   studlab = studylab,
                   method="Inverse",
                   byvar = filtered$organ,
                   comb.random = T,
                   comb.fixed = F,
                   title=x)
  
  png(sprintf("/home/chit/Desktop/Hiwi_BMC/Overlapstudies/%s_nonocular.png", x), width = 1000, height = 1200)
  forest(meta)
  #grid.text(x, .5, 1-1/sizefactor, gp=gpar(cex=2))
  dev.off()
}



