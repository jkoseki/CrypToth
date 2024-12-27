### Estimate average scores
##
#
library(tidyr)
library(data.table)
library(stringr)


f.list     <-  list.files("./", full.names = TRUE, recursive = TRUE)
change.id  <-  grep("-Spot_Scores.csv", f.list)
f.list     <-  f.list[change.id]

# Set your target PDB IDs
pdb.id     sapply(strsplit(f.list, ".//"), function(x){x[2]})
pdb.id     sapply(strsplit(pdb.id, "/"),   function(x){x[1]})

# The six Probes used within msmd are used in the following notation.
## A00 :  Benzene
## A01 :  Isopropanol
## A20 :  Phenol
## A37 :  Imidazole
## B71 :  Acetonitrile
## E20 :  Ethylene glycol
sol.id     <-  c('A00', 'A01', 'A20', 'A37', 'B71', 'E20')

for (PN in pdb.id) {
  tf.list  <-  f.list[grep(PN, f.list, ignore.case=TRUE)]
  
  for (SLV in sol.id) {
    tf2.list  <-  tf.list[grep(SLV, tf.list)]
    
    for (fn in 1:length(tf2.list)) {
      if (fn != 1) {
        RS         <-  read.csv(tf2.list[fn])
        AS[, 3:5]  <-  AS[, 3:5] + RS[, 3:5]
      } else {
        AS <- read.csv(tf2.list[fn])
      }
    }
    
    AS[, 3:5]  <-  AS[, 3:5] / length(tf2.list)
    
    out.name   <-  paste0(PN, "_", SLV, "-Average-Scores.csv")
    write.csv(AS, file = out.name, row.names = FALSE, quote = FALSE)
  }
}




