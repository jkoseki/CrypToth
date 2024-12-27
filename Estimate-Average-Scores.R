### Estimate average scores
##
#
library(tidyr)
library(data.table)
library(stringr)


f.list     <-  list.files("./", full.names = TRUE, recursive = TRUE)
change.id  <-  grep("-Spot_Scores.csv", f.list)
f.list     <-  f.list[change.id]

Of.list    <-  list.files("../", full.names = TRUE, recursive = TRUE)
change.id  <-  grep("-Spot_Scores.csv", Of.list)
Of.list    <-  Of.list[change.id]

change.id  <-  grep("Additional-TDA", Of.list)
Of.list    <-  Of.list[-change.id]

change.id  <-  grep("old-data", Of.list)
Of.list    <-  Of.list[-change.id]
change.id  <-  grep("4P0I", Of.list)
Of.list    <-  Of.list[-change.id]


f.list     <-  c(Of.list, f.list)


pdb.id     <-  c('1FVR', '1FXX', '1NEP', '2AM9', '2W9T', '3NX1', '3P53', '3QXW')
#pdb.id     <-  c('1JWP')
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




