#########
## Genetic distance
########
##  Genetic distance
########

# sq = DNAStringSet
# samples_id = samples_id or names in the id column to select the haplotypes beloning to a given sample.
# db = dataframe that contains at least 3 columns, id, hap, read_num.

genetic_dist = function(samples_id, db, sq){
  L = width(sq)[1]
  ## Consensus matrix per sample.
  cm_id = list()
  for(i in samples_id){
    temp = matrix(0, 4, L)
    db_set = db[db$id == i, c('hap','read_num')]
    for(h in 1:nrow(db_set)) temp = temp + consensusMatrix(sq[db_set$hap[h]])[1:4,] * db_set$read_num[h]
    cm_id[[i]] = temp / sum(temp[,1])
  }
  ## Compute distance
  d = data.frame()
  for(a in 1:length(cm_id)){
    for(b in 1:a){
      temp = data.frame(x=samples_id[a], y=samples_id[b],
                        g_dist = sqrt( sum( sapply(1:L, function(i) sum((cm_id[[a]][,i] - cm_id[[b]][,i])^2) ) )/L )     )
      d = rbind(d, temp)
    }
  }
  return(d)
}
