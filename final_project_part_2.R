#II-1 what is the family of the species in your data set? describe some basic
#biological information about the family
#misumenops is a genus of spiders that are in the family Thomisidae, or 
#crab spiders. Crab spiders are ambush predators and do not spin webs like other 
#spider species. The front two legs of spiders in this family are longer and
#larger than their back legs, and their bodies are flat, making them resemble 
#crabs. 
#https://www.inaturalist.org/taxa/47866-Thomisidae

#II-2 what are the uses of a sequence alignment? A sequence alignment is used to
#determine how similiar closely related species are, and can also show 
#evolutionary lines.

#II-3 find one or more scientific papers that have used a sequence alignment of 
#DNA, RNA, or amino acids to perform some task 
#Wang, Y., Wu, H. & Cai, Y. A benchmark study of sequence alignment methods for protein clustering. 
#BMC Bioinformatics 19 (Suppl 19), 529 (2018). https://doi.org/10.1186/s12859-018-2524-4
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2524-4#citeas

library("rentrez")
library("stringr")
library("ape")
library(Biostrings)
library("zoo")
mis_seq <- read.csv(file = "Misumenops_R16S-ND1.csv")
#II04 there are 14 species in the data set 
nrow(mis_seq)
#II-5 download all the sequences to the same fasta file 
seq_fasta <- vector(mode = "character", length = "14")
for (it in 1:14) {
  print(it)
  accn_search <- entrez_search(db="nucleotide", term = mis_seq[it,2])
  seq_fasta[it] <- entrez_fetch(db="nucleotide", id = accn_search$ids, 
                                rettype = "fasta")
}
write(x=seq_fasta, file = "misumenops_all_sequences.fasta")
mis <- readDNAStringSet(filepath = "misumenops_all_sequences.fasta")
#II-7 Use the mafft website to preform sequence alignment 
#II-8 rename each sequence in the alignment 
mis_align <- readDNAStringSet(filepath = "misumenops alignment.fasta")
names(mis_seq)
names(mis_align)
names(mis_align) <- mis_seq[,1]
names(mis_align)
#II-9 Compare the length of each sequence before and after alignment
mis_mat <- as.matrix(mis_align)
length(mis_align)
nchar(mis_align)
nchar(mis)
#each sequence varied in length before the alignment, but after the alignment
#the length of each sequence is 585. The shortest sequence before alignment was
#576 nucleotides long and the longest was 583.
#II-10 the species that showed greatest increase in sequence length 
for (it in 1:14){
  if(nchar(mis[it])==576){
    print(mis_seq[it,1])
  }
}
#II-11 the sequences that have the largest number of gaps
alpha_freq <- alphabetFrequency(mis_align, as.prob = FALSE)
alpha_freq
gaps <- alpha_freq[,"-"]
max(gaps)
gaps == max(gaps)
seq_names <- names(mis_align)
seq_names[gaps==max(gaps)]
#II-12 GC content of each sequence 
help("GC.content")
GC <- vector(mode = "numeric", length = 14)
for (it in 1:14) {
  GC[it] <- (alpha_freq[it,"G"]+alpha_freq[it,"C"]/585)
  print(GC[it])
}
#II-13 histogram of GC content across each sequence
hist(GC, xlab = "GC content of each sequence")
#the plot shows one extreme value of 51%, with the majority of 
#sequences having a GC content between 52%-56%
#II-14 compute proportion of gaps in each sequence 
gap_prop <- vector(mode = "numeric",length = 14)
for (it in 1:14) {
  gap_prop[it] <- (alpha_freq[it,"-"]/585)
  print(gap_prop[it])
}
#II-15 Produce the consenus matrix 
con_mat <- consensusMatrix(x=mis_align)
#II-16 Make a plot of distribution of proportion of gaps.
#the gaps are mostly concentrated in the beginning of the sequence around 
#position 100-200
con_gap <- rollmean(x=con_mat["-",], k=10)
plot(con_gap, xlab = "Rolling mean of gaps per loci")
#II-17 Find alignment positions without any gaps 
sum(con_mat["-",] == 0)
sum(con_mat)
sum(con_mat) - sum(con_mat["-",])

#II-19 Construct random DNA sequence wiwth length equal to length of longest sequence
nuc <- c("A","T","G","C")
random_DNA <- sample(nuc,size = 585, replace = TRUE)
paste(random_DNA, collapse = "")
fake_alignment <- readDNAStringSet(filepath = "fake_dna_alignment.fasta")
length(fake_alignment)   
#II-22 The length of the sequences in the new alignment is 
width(fake_alignment)  
#the length of the sequences in the old alignment was 
width(mis_align)
#II-23 the inclusion of the random sequence changed the alignment because the 
#program was trying to force similarities between the real species and the fake 
#dna that would have no evolutionary similarities to the previous data, so it 
#would put in more gaps to compensate 
#II-24 Produce a consensus matrix for new data
fake_consensus <- consensusMatrix(x=fake_alignment)
#plot distribution of proprtion of gaps
fake_gap <- rollmean(x=fake_consensus["-",], k=10)
plot(fake_gap, xlab = "Rolling mean of gaps per loci")
#the x-axis of this plot is showing the positions on the dna from 0 to 772. 
#the y axis is showing the average number of gaps at each position. 
#since we added fake data to the alignment, there is a large number of gaps for 
#the first 0-350 positions in the dna, as over half of the sequence has some number
#of gaps. 
                                   
                                   
                                   
                                   
                                  