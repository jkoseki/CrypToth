### Searching common atoms
##
#
library(tidyr)
library(data.table)
library(stringr)

f.list     <-  list.files("./", full.names = TRUE, recursive = TRUE)
org.list   <-  f.list
change.id  <-  grep("-Change-point.csv", f.list)
f.list     <-  f.list[change.id]

d.list     <-  sapply(strsplit(f.list, ".//"), function(x){x[2]})
case.list  <-  sapply(strsplit(d.list, "/"), function(x){x[1]})
solv.list  <-  sapply(strsplit(case.list, "_"), function(x){x[2]})


for (i in 1:length(case.list)) {
  f.solv.list  <-  f.list[i]
  d.name       <-  paste0(case.list[i], "/Pickup-Variation-Data")
  d.solv       <-  sapply(strsplit(d.name, "/"), function(x){x[1]})
  sc.list      <-  org.list
  score.id     <-  grep(d.name, sc.list)
  sc.list      <-  sc.list[score.id]

  sc.list  <-  sapply(strsplit(sc.list, ".pdf"), function(x){x[1]})
  sc.list  <-  sapply(strsplit(sc.list, "Pickup-Variation-Data/"), function(x){x[2]})  
  sc.list  <-  sapply(strsplit(sc.list, "WT-"), function(x){x[2]})
  sc.list  <-  sapply(strsplit(sc.list, "_"  ), function(x){x[2]})
  sc.list  <-  sc.list[-length(sc.list)]

  sc.id    <-  which(abs(as.numeric(sc.list)) >= 5)
  
  TDA.point     <-  read.csv(f.solv.list)
  PU.TDA.point  <-  TDA.point[sc.id, c(3, 4, 6)]
  PU.TDA.point  <-  cbind(PU.TDA.point, as.vector(as.numeric(sc.list)[sc.id]))
  names(PU.TDA.point)[4]  <-  "ID Score"

  out.name  <-  paste0(d.solv, "-CP.csv")
  write.csv(PU.TDA.point, file = out.name, quote = FALSE, row.names = FALSE) 
}


s.list     <-  list.files("./spots/", full.names = TRUE)
s.list.id  <-  grep("U-spots", s.list)
s.list     <-  s.list[s.list.id]


cp.list    <-  list.files("./", full.names = TRUE)
cp.list.id <-  grep("-CP.csv", cp.list)
cp.list    <-  cp.list[cp.list.id]

num.cp     <-  length(cp.list)

for (i in 1:num.cp) {
  cp.data  <-  read.csv(cp.list[i])

  if (nrow(cp.data) != 0) {
    r.num    <-  nrow(cp.data)
    
    cp.data.at  <-  strsplit(cp.data$PDB.Number, split = ";")
    cp.data.sc  <-  cp.data.at
    
    for (j in 1:r.num) {
      cp.data.sc[[j]][1:length(cp.data.at[[j]])]  <-  cp.data$ID.Score[j]
    }
    
    sc.sum    <-  cbind(as.numeric(unlist(cp.data.at)), as.numeric(unlist(cp.data.sc)))
    sc.sum    <-  sc.sum[order(sc.sum[, 1]), ]
    
    sol.name  <-  sapply(strsplit(cp.list[i], "-CP.csv"), function(x){x[1]})
    sol.name  <-  sapply(strsplit(sol.name,   "//"),      function(x){x[2]})  
    
    colnames(sc.sum)  <-  c(colnames(cp.data)[3:4]) 
    
    out.name  <-  paste0(sol.name, "-Score.csv")
    write.csv(sc.sum, file = out.name, quote = FALSE, row.names = FALSE) 
  }
}


sol.sc.list  <-   sapply(strsplit(cp.list, "-CP.csv"), function(x){x[1]})
sol.sc.list  <-   sapply(strsplit(sol.sc.list, "//" ), function(x){x[2]})

num.spot   <-  length(s.list)


score.list <-  list.files("./")
id.sc.list <-  grep("-Score.csv", score.list)
score.list <-  score.list[id.sc.list]



for (i in 1:num.spot) {
  spot.dat  <-  read.table(s.list[i])
  count.and.score  <-  data.table(matrix(data = 0, nrow = nrow(spot.dat), ncol = 2))
  spot.dat  <-  cbind(spot.dat, count.and.score)
  colnames(spot.dat)  <-  c("Res.name", "PDB.Number", "Count", "Ave.score")
  
  org.dat   <-  spot.dat
  
  for (j in 1:num.cp) {
    in.name  <-  paste0(sol.sc.list[j], "-Score.csv")

    if (length(which(score.list == in.name)) != 0) {
      t.dat    <-  read.csv(in.name)
      
      spot.dat <- org.dat
      
      for (k in 1:nrow(t.dat)) {
        spot.id  <-  which(spot.dat$PDB.Number == t.dat$PDB.Number[k])
        
        if (length(spot.id) > 0) {
          spot.dat[spot.id, ]$Ave.score  <-  (spot.dat[spot.id, ]$Ave.score * spot.dat[spot.id, ]$Count)
          spot.dat[spot.id, ]$Ave.score  <-   spot.dat[spot.id, ]$Ave.score + abs(t.dat[k, 2])
          spot.dat[spot.id, ]$Count      <-   spot.dat[spot.id, ]$Count     +  1
          spot.dat[spot.id, ]$Ave.score  <-   spot.dat[spot.id, ]$Ave.score / spot.dat[spot.id, ]$Count
        }
      }
      
      out.name  <-  paste0("Spot", str_pad(i, width = 2, pad = "0"), "-", sol.sc.list[j], ".csv")
      write.csv(spot.dat, file = out.name, quote = FALSE, row.names = FALSE) 
    }
  }
}


for (i in 1:num.cp) {
  t.sp.list  <-  list.files("./", full.names = TRUE)
  
  t.id.sol   <-  grep(sol.sc.list[i], t.sp.list)
  t.id.spot  <-  grep("Spot",         t.sp.list)
  t.sp.list  <-  t.sp.list[intersect(t.id.sol, t.id.spot)]
  t.id.spot  <-  grep("Scores",       t.sp.list)
  
  if (length(t.id.spot) > 0) {
    t.sp.list  <-  t.sp.list[-t.id.spot]
  }

  if (length(t.sp.list) > 0) {
    Spot.name  <-  sapply(strsplit(t.sp.list,  "-"),  function(x){x[1]})
    Spot.name  <-  sapply(strsplit(Spot.name,  "//"), function(x){x[2]})
    
    Spot.score <-  data.table("Spot"      = Spot.name,
                              "Res.num"   = 0,
                              "Tot.Count" = 0,
                              "Ave.Score" = 0,
                              "Cor.Score" = 0)
    
    for (j in 1:length(t.sp.list)) {
      spot.dat  <-  read.csv(t.sp.list[j])
      
      Spot.score[j, 2]  <-  nrow(spot.dat)
      Spot.score[j, 3]  <-  sum(spot.dat$Count)
      if (Spot.score[j, 3] > 0) {
        Spot.score[j, 4]  <-  sum(spot.dat$Count * spot.dat$Ave.score) / Spot.score[j, 3]
        Spot.score[j, 5]  <-  Spot.score[j, 3]   * Spot.score[j, 4]    / Spot.score[j, 2]
      }
    }
    
    out.name  <-  paste0(sol.sc.list[i], "-Spot_Scores", ".csv")
    write.csv(Spot.score, file = out.name, quote = FALSE, row.names = FALSE)
  }
}









