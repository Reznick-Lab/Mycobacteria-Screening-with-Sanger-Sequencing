###--TRANSFORMING CSV DATASET INTO CATEGORICAL VARIABLES--##

##downloading libraries and packages
library(dplyr)

##set working directory into the file you want accessed setwd()

##replacing using a lookup table

# 1. Your lookup table (can also be read from a CSV)

df<- read.csv("combined_blast_results.csv")

lookup <- data.frame(
  subject_id = c(
    "AM396443.1", 
    "AM396444.1", 
    "AM396445.1", 
    "AM396447.1",
    "AM396448.1",
    "AM396449.1",
    "AM902926.1",
    "AM902929.1",
    "AM902935.1"
    ),
  species = c(
    "Mycobacterium chelonae ITS1, isolate 29/02",
    "Mycobacterium chelonae ITS1, isolate T5-2",
    "Mycobacterium chelonae ITS1, isolate T14",
    "Mycobacterium marinum ITS1, isolate 2B",
    "Mycobacterium marinum ITS1, isolate 9/04",
    "Mycobacterium marinum ITS1, isolate 131/03",
    "Mycobacterium fortuitum ITS1, strain S7",
    "Mycobacterium fortuitum ITS1, strain 277/2/01",
    "Mycobacterium fortuitum ITS1, strain 276/3/01"
  )
)


# 2. Merge and replace the column
df <- df %>%
  left_join(lookup, by = "subject_id") %>%
  mutate(subject_id = ifelse(is.na(species), subject_id, species)) %>%  # Replace if found
  select(-species)  # Remove helper column

##use head(df) to inspect joined data frame


##write new file
write.csv(df, "species_id_blast_results.csv", row.names =FALSE )

