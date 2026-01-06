#SANGER SEQUENCING THROUGH R--DESIGNED FOR MYCOBACTERIA ITS SEQUENCES#

##first set up all your packages##
install.packages("BiocManager")
library(BiocManager)
BiocManager::install(c("sangerseqR", "annotate", "Biostrings", "rBLAST", "bio3d", "seqinr")) 
##packages associated with Biocmanager should be downloading at this step. 
library(sangerseqR)
##create a directory (folder) that accesses your .ab1 files
##getwd() ##pull up your working directory, what file are you in? 
library(annotate)##command for accessing annotate function from the Bioconductor package
##this accesses tools for conducting a BLAST search in R 
library(Biostrings)
##accessing rBLAST library to use tools like BLASTn
library(rBLAST)
library(bio3d)
library(seqinr)
getwd()##figure out which file you are in

##--now that you have your appropriate libraries ready, you can start working with your sequences--##
##setwd()-> set working directory into file: "4-11-25-10_MYCOITS-VP-DS"
ITS<-read.abif("F11L_A01_015.ab1") ##1st well on the plate, this function reads the ab1 file
##then we convert the ab1 file into a Sangerseq object (from the library sangerseqR package you downloaded earlier)

##4-11-25-10_MYCOITS-VP-DS

#----LOOP THROUGH ALL FILES AND CONVERT TO FASTA----#
##set working directory to .ab1 files
setwd("/Users/froglord/Desktop/4-11-25-10_MYCOITS-VP-DS")

##Create a folder for your fasta files to go to
output_directory<-"converted_fasta2"
if(!dir.exists(output_directory)) {
  dir.create(output_directory)
}

##get list of ab1 files
ab1_files <- list.files(pattern = "\\.ab1$")

##for loop to process all files and convert into a FASTA format
for (file in ab1_files) { 
  ##read chromatogram 
  abif_data<-read.abif(file)
  ##extract sequence
  sequence<- sangerseq(abif_data)
  #get primary DNA sequence
  dna_seq<- primarySeq(sequence)
  #convert to DNAStringSet -> fasta
  fasta<- DNAStringSet(dna_seq)
  #name sequence without .ab1 file
  sample_name <- tools::file_path_sans_ext(basename(file))
  names(fasta)<- sample_name
  
  #write to .fasta file
  fasta_file<- paste0(sample_name,".fasta")
  
  #save into output folder
  fasta_path<- file.path(output_directory, fasta_file)
  writeXStringSet(fasta, filepath = fasta_path)
  ##prints status message that keeps track of what has been converted
  cat("✅ Converted:", file, "→", fasta_path, "\n")
}

##--find all your fasta files in converted_fasta file!--##

#------TO SEND A QUERY TO BLAST THROUGH BIOCONDUCTOR-------##


##--first make a local database below--##
# 1. Install needed packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("rentrez")
library(rentrez)

# 2. Define the accession numbers you want to blast against
accessions <- c(
  "AM396443.1", #Mchelonae ITS1, 29/02
  "AM396444.1", #Mchelonae ITS1, T5-2
  "AM396445.1", #Mchelonae ITS1, T14
  "AM396447.1", #Mmarinum ITS1, 2B
  "AM396448.1", #Mmarinum ITS1, 9/04
  "AM396449.1", #Mmarinum ITS1, 131/03
  "AM902926.1", #Mfortuitum ITS1, S7
  "AM902929.1", #Mfortuitum ITS1, 277/2/01
  "AM902935.1"  #Mfortuitum ITS1, 276/3/01
)  

##set working directory to coverted_fasta folder 

# 3. Fetch sequences and write them to one FASTA file
output_fasta <- "reference_sequences.fasta"
all_seqs <- sapply(accessions, function(acc) {
  cat("Downloading:", acc, "\n")
  entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
})
writeLines(unlist(all_seqs), output_fasta)
cat("All sequences written to", output_fasta, "\n")

#build a local BLAST database
db_name <- "ITS_custom_db"
system(paste("makeblastdb -in", output_fasta, "-dbtype nucl -out", db_name))
cat("BLAST database created:", db_name, "\n")

# Set path to your reference FASTA file
ref_fasta <- "reference_sequences.fasta"

# Run makeblastdb via system call
make_db_cmd <- paste("makeblastdb -in", ref_fasta, "-dbtype nucl -out", db_name)
system(make_db_cmd)

##--create a loop to query multiple fasta files--##

#set paths
query_dir <- "queries/" ##this is where all your queried Fasta files will go
db_name <- "ITS_custom_db" 
blast_output <- "blast_results/" #where your blast results will go

dir.create("blast_results", showWarnings = FALSE)


##loop and blast each query file 
for (query_path in fasta_files) {
  query_name <- tools::file_path_sans_ext(basename(query_path))
  output_path <- file.path(blast_output, paste0(query_name, "_blast.txt"))
  
  blast_cmd <- paste("blastn",
                     "-query", shQuote(query_path),
                     "-db", shQuote(db_name),
                     "-out", shQuote(output_path),
                     "-outfmt 6")  # tabular format (can customize)
  
  cat("🔍 Running BLAST for:", query_name, "\n")
  system(blast_cmd)
  cat("✅ Results saved to:", output_path, "\n\n")
}

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# Set path to BLAST results folder
blast_folder <- "blast_results/"
blast_files <- list.files(blast_folder, pattern = "_blast\\.txt$", full.names = TRUE)

# Function to read one file and add query name
read_blast <- function(file) {
  df <- read.delim(file, col_names = FALSE, show_col_types = FALSE)
  colnames(df) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  df$query_file <- tools::file_path_sans_ext(basename(file))
  return(df)
}

# Combine all BLAST results
blast_data <- bind_rows(lapply(blast_files, read_blast))


files <- list.files("blast_results", pattern = "\\.txt$", full.names = TRUE)

for (file in files) {
  blast_data <- read_tsv(file, col_names = FALSE, show_col_types = FALSE)
  colnames(blast_data) <- c(
    "query_id", "subject_id", "percent_identity", "alignment_length",
    "mismatches", "gap_opens", "q_start", "q_end",
    "s_start", "s_end", "evalue", "bit_score"
  )
  
  out_file <- sub("\\.txt$", ".csv", file)
  write.csv(blast_data, out_file, row.names = FALSE)
  cat("✅ Saved:", out_file, "\n")
}

species <- c("Mycobacterium chelonae ITS1 isolate 29/02",
             "Mycobacterium chelonae ITS1 isolate T5-2",
             "Mycobacterium chelonae ITS1 isolate T14",
             "Mycobacterium marinum ITS1 isolate 2B",
             "Mycobacterium marinum ITS isolate 9/04",
             "Mycobacterium marinum ITS1 isolate 131/03",
             "Mycobacterium fortuitum ITS1 strain S7",
             "Mycobacterium fortuitum ITS1 strain 277/2/01",
             "Mycobacterium fortuitum ITS1 strain 276/3/01")

species_lookup <- data.frame(
  subject_id = "accessions",
  species = "species"
)
##organizing individual blast results into one big txt file
blast_txt_files <- list.files("blast_results", pattern = "\\.txt$", full.names = TRUE)

blast_result <- read_tsv("blast_txt_files", col_names = FALSE, show_col_types = FALSE)

blast_all <- lapply(blast_txt_files, function(file) {
  data <- read_tsv(blast_txt_files, col_names = FALSE, show_col_types = FALSE)
  data$file <- basename(file)  # Add filename for traceability
  return(data)
}) |> bind_rows()

##add column names
colnames(blast_all)[1:12] <- c(
  "query_id", "subject_id", "percent_identity", "alignment_length",
  "mismatches", "gap_opens", "q_start", "q_end",
  "s_start", "s_end", "evalue", "bit_score"
)

setwd("/Users/froglord/Desktop/SANGER-SEQUENCING/SANGER_TEST/converted_fasta")
blast_txt_files <- list.files("blast_results", pattern = "\\.txt$", full.names = TRUE)
