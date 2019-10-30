library(dplyr)
library(tidyr)
library(ggplot2)
library(metafor)
library(meta)
library(dmetar)
library(grid)

classi <- /path/to/table/of/classification
overlap <- /path/to/raw/data

#adjust the scale of the data to patients
overlap <- overlap %>%
  mutate(total=ifelse(scale=="eyes", total/2, total),
         case=ifelse(scale=="eyes", case/2, case)) %>%
  mutate(case=ifelse(scale=="percentage", total*(case/100), case)) %>%
  mutate(studylab=paste0(study,"_",autoimmune))

##meta-analysis for each autoimmune disease
###save the result value for each input
seTE <- c() #standard error
TE <- c() #proportion
lower <- c() #lower confidence interval
upper <- c() #upper confidence interval
validauto <- c() #autoimmune disease name
classification <- c() #classification of it

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
  
  png(/path/to/save/the/image, width = 1000, height = 1200)
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

png(/path/to/save/the/image, width = 800, height = 700)
forest(meta_autoall,
       comb.random = T,
       comb.fixed = F,
       leftcols = "studlab",
       colgap.forest.left = "40mm",
       col.diamond.random = "lightblue",
       col.square = "darkblue",
       col.square.lines = "darkblue")
dev.off() ##saved


#meta-analysis for non-ocular manifestation
nonocular <- overlap %>%
  filter(ocular=="no")

##same as before
nseTE <- c()
nTE <- c()
nlower <- c()
nupper <- c()
nvalidauto <- c()
nclassification <- c()
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
  
  png(/path/to/save/the/image, width = 1000, height = 1200)
  forest(meta)
  #grid.text(x, .5, 1-1/sizefactor, gp=gpar(cex=2))
  dev.off() ##saved
  
  nseTE <- c(nseTE, meta$seTE.random)
  nTE <- c(nTE, meta$TE.random)
  nlower <- c(nlower, meta$lower.random)
  nupper <- c(nupper, meta$upper.random)
  nvalidauto <- c(nvalidauto, x)
  
  y <- classi[classi$ad==x,]$Classification
  nclassification <- c(nclassification, y)
}

nautoall <- data.frame("seTE"=meta:::backtransf(nseTE, sm="PLOGIT"), 
                       "TE"=meta:::backtransf(nTE, sm="PLOGIT"), 
                       "lower"=meta:::backtransf(nlower, sm="PLOGIT"),
                       "upper"=meta:::backtransf(nupper, sm="PLOGIT"),
                       nvalidauto, "classification"=nclassification)



meta_nautoall <- metagen(TE=TE,
                         seTE=seTE,
                         lower = lower,
                         upper = upper,
                         data=nautoall,
                         studlab = nvalidauto,
                         byvar = nautoall$classification,
                         sm="Prevalence")

png(/path/to/save/the/image, width = 800, height = 700)
forest(meta_nautoall,
       comb.random = T,
       comb.fixed = F,
       leftcols = "studlab",
       colgap.forest.left = "40mm",
       col.diamond.random = "lightblue",
       col.square = "darkblue",
       col.square.lines = "darkblue")
dev.off() ##saved

