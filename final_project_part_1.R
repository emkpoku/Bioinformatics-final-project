#I-1 an operon is a group of genes that share one promoter 
#I-2 the e.coli lac operon sequence encodes for three genes
#I-3 the main function of the lac operon sequence is to. allow the bacteria to 
#break down lactose so that it can be used as energy
# source: https://www.khanacademy.org/science/ap-biology/gene-expression-and-regulation/regulation-of-gene-expression-and-cell-specialization/a/the-lac-operon
install.packages("rentrez")
library("rentrez")
install.packages("stringr")
library("stringr")
library(Biostrings)
#I-4 download the DNA sequence for E.coli operon 
acc_code <- "J01636.1"
accn_search <- entrez_search(db="nucleotide", term = acc_code)
accn_search
seq_fasta <- entrez_fetch(db="nucleotide", id = accn_search$ids, 
        rettype = "fasta")
seq_fasta
#I-5 save sequence as FASTA file 
write(x = seq_fasta, file ="ecoli_J01636_1.fasta" )
#change the sequence into a string of characters
lac_operon <- as.character(seq_fasta)[[1]]
class(lac_operon)
str_length(lac_operon)
#I-6 isolate the gene sequences for lacI, lacZ,lacY, and lacA
lacI_rep <- substr(x = lac_operon, start = 79, stop = 1161)
lacZ <- substr(x = lac_operon, start = 1284, stop = 4358)
lacY <- substr(x = lac_operon, start = 4410, stop = 5663)
lacA <- substr(x = lac_operon, start = 5727, stop = 6338)
#I-7 isolate the introns in the gene 
intron1 <- substr(x = lac_operon, start = 1162, stop = 1283)
intron2 <- substr(x = lac_operon, start = 4359, stop = 4409)
intron3 <- substr(x = lac_operon, start = 5664, stop = 5726)
#I-8 compute the GC content for each of the genes 
lacIgc <- str_count(string = lacI_rep, pattern = "G") + 
  str_count(string = lacI_rep, pattern = "C")
#the GC content of lacI is 
lacIgc <- lacIgc/str_length(lacI_rep)

lacZgc <- str_count(string = lacZ, pattern = "G") + 
  str_count(string = lacZ, pattern = "C")
#the GC content of lacZ is 
lacZgc <- lacZgc/str_length(lacZ)

lacYgc <- str_count(string = lacY, pattern = "G") + 
  str_count(string = lacY, pattern = "C")
#The GC content of lacY is 
lacYgc <- lacYgc/str_length(lacY)
  
lacAgc <- str_count(string = lacA, pattern = "G") + 
  str_count(string = lacA, pattern = "C")
#the GC content of lacA is 
lacAgc <- lacAgc/str_length(lacA)

#I-9 compute the gc content for each of the intron sequences 
intron1gc <- str_count(string = intron1, pattern = "G") + 
  str_count(string = intron1, pattern = "C")
#the GC content of intron 1 is 
intron1gc <- intron1gc/str_length(intron1)

intron2gc <- str_count(string = intron2, pattern = "G") + 
  str_count(string = intron2, pattern = "C")
#the GC content of intron 2 is 
intron2gc <- intron2gc/str_length(intron2)

intron3gc <- str_count(string = intron3, pattern = "G") + 
  str_count(string = intron3, pattern = "C")
#the GC content of intron 3 is 
intron3gc <- intron3gc/str_length(intron3)

#I-10 which sequence region (genes or introns) has the highest GC content 
intronGC <- intron1gc+intron2gc+intron3gc
geneGC <- lacIgc+lacZgc+lacYgc+lacAgc
intronGC>geneGC
geneGC>intronGC
#the genes have the highest GC content 
str_length(lacI_rep)
str_length(lacZ)
str_length(lacY)
str_length(lacA)
#I-11 the lacZ gene has the longest gene sequence
#I-12 create a new sequence by binding together the coding regions 
lac_operon_CDS <- c(lacZ,lacY,lacA)
lac_operon_CDS <- paste(lac_operon_CDS, collapse = "")
lac_operon_CDS

#I-13 create a new sequence that binds together only the non coding regions 
lac_non_coding <- c(intron1, intron2, intron3)
lac_non_coding <- paste(lac_non_coding, collapse = "")

#I-14 find the researchers short read in a gene of the lac operon 
short_read <- "GTCGGCATCATGTTCACCAT"
str_count(string=lacI_rep, pattern = short_read)
str_count(string = lacZ, pattern = short_read)
#the researchers short read was found in gene lacY 
str_count(string=lacY, pattern = short_read)
str_count(string = lacA, pattern = short_read)

#I-15 compute the number of each type of base pair necessary 
full_lacZ <- substr(x = lac_operon, start = 1246, stop = 4358)
#number of A's
countA <- str_count(string = full_lacZ, pattern = "A")
#number of U'S
countU <- str_count(string = full_lacZ, pattern = "T")
#number of C's
countC <- str_count(string = full_lacZ, pattern = "C")
#number of G's
countG <- str_count(string = full_lacZ, pattern = "G")

#I-16 count the number of each amino acid produced by the sequence of the 
#lacZ gene 
lacZ_dna <- readDNAStringSet(filepath = "ecoli_J01636_1.fasta")
lacZ_dna <- DNAStringSet(x=lacZ_dna,start = 1284, end = 4358)
lacZ_AA <- translate(lacZ_dna)
alphabetFrequency(lacZ_AA)





