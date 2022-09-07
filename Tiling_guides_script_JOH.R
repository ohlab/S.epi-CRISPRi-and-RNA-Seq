############### STRAND 5'-3' ############

# load packages 
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))

# change working directory 
setwd ('/Users/')


### FUNCTIONS ####


#function to find PAM site (top strand)
FindCC<-function (Sequences) 
{str_locate_all(Sequences,"CC")}


#reverse comp guide 
ReverseCompGuide<-function(x) {
  str_c(rev(comp(s2c(x))), collapse = "")
}


#remove NGG/PAM sequence
RemoveNGG<-function(x) { 
  substr(x,1,(nchar(x) - 3))}


###############################

# Input the genome sequence 

genome <- readDNAStringSet(filepath = "Tu3298_whole_genome_singleseq.fasta")


# Find PAM sites (searching for 'CC') across the genome (bottom strand)
genome_rev_comp <- reverseComplement(genome) # reverse complement genome 
GuideLocations_rev <- lapply(genome_rev_comp, FindCC)
head(GuideLocations_rev)
GuideLocationsdf_rev <- data.frame(GuideLocations_rev)
colnames(GuideLocationsdf_rev) <- c("start", "end")

# Making guides from "CC" sites 
MakeGuidesList_rev <- function(r,N) {
  Start <- as.numeric(r[["start"]])
  Guide <- substr(genome_rev_comp, Start, Start+N)
}

OPTS.N <- 20 #guide length
AllGuidesList_rev <- apply( GuideLocationsdf_rev, 1, MakeGuidesList_rev, OPTS.N) # applying function to make guides of length 20 


# write out 
write.csv(AllGuidesList_rev, file = "Tiling_bottomstrand_ALL.GUIDES.csv")


# Reverse-comp guides (b/c made from "CC" site on top strand; convert to guide upstream of "GG site on bottom strand)
CorrectedGuidesList_rev <- lapply(AllGuidesList_rev, ReverseCompGuide)
# setting to full guides list b/c that will be the name once the bottom strand guides are appended 
Full_guides_list <- CorrectedGuidesList_rev

# Calculating GC content 
Guides_as_char <- lapply (Full_guides_list, s2c)
With_GC_content <- data.frame(cbind(lapply(Guides_as_char, GC),
                                    Full_guides_list))

colnames(With_GC_content) <- c("GC_perc", "Guide")
With_GC_content$Guide_list <- lapply(With_GC_content$Guide, RemoveNGG)

# subset out guides with GC content greater than 35 % 
Final_guides <- subset(With_GC_content, GC_perc > .35, select=c(Guide_list))
Final_guides_matrix <- as.matrix(Final_guides)
write.csv(Final_guides_matrix, file = "FINAL.GUIDES_bottomstrand.csv")



############### STRAND 3'-5' #############

# load packages 
suppressPackageStartupMessages(library(BiocInstaller))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(refGenome))
suppressPackageStartupMessages(library(genbankr))
suppressPackageStartupMessages(library(rentrez))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(annotate))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringi))

setwd ('/Users/spotom/Desktop/Tiled_guides')


### FUNCTIONS ####


#function to find PAM site (top strand)
FindCC<-function (Sequences) 
{str_locate_all(Sequences,"CC")}


#reverse comp guide 
ReverseCompGuide<-function(x) {
  str_c(rev(comp(s2c(x))), collapse = "")
}


#remove NGG/PAM sequence
RemoveNGG<-function(x) { 
  substr(x,1,(nchar(x) - 3))}

#Remove_GC<-function(x) { 
#if (x > .35) { 
#print(x)}
#else {
#print(NA)}
#}

###############################

# Input the genome sequence 

genome <- readDNAStringSet(filepath = "Tu3298_whole_genome_singleseq.fasta")

# Find PAM sites (searching for 'CC') across the genome (top strand)
GuideLocations <- lapply(genome, FindCC)
head(GuideLocations)
GuideLocationsdf <- data.frame(GuideLocations)
colnames(GuideLocationsdf) <- c("start", "end")

# Making guides from "CC" sites 
MakeGuidesList <- function(r,N) {
  Start <- as.numeric(r[["start"]])
  Guide <- substr(genome, Start, Start+N)
}

OPTS.N <- 20 # guide length
AllGuidesList <- apply(GuideLocationsdf, 1, MakeGuidesList, OPTS.N) # applying function to make guides of length 20 
df_AllGuidesList <- data.frame(AllGuidesList, stringsAsFactors = FALSE)

# write out 
write.csv(AllGuidesList, file = "Tiling_topstrand_ALL.GUIDES.csv")


# Reverse-comp guides (b/c made from "CC" site on top strand; convert to guide upstream of "GG site on bottom strand)
CorrectedGuidesList <- lapply(AllGuidesList, ReverseCompGuide)
# setting to full guides list b/c that will be the name once the bottom strand guides are appended 
Full_guides_list <- CorrectedGuidesList

# Calculating GC content 
Guides_as_char <- lapply (Full_guides_list, s2c)
With_GC_content <- data.frame(cbind(lapply(Guides_as_char, GC),
                                    Full_guides_list))

colnames(With_GC_content) <- c("GC_perc", "Guide")
With_GC_content$Guide_list <- lapply(With_GC_content$Guide, RemoveNGG)

# subset out guides with GC content greater than 35 % 
Final_guides <- subset(With_GC_content, GC_perc > .35, select=c(Guide_list))
Final_guides_matrix <- as.matrix(Final_guides)
write.csv(Final_guides_matrix, file = "FINAL.GUIDES_topstrand.csv")


