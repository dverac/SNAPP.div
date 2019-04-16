##############################
## FUNCTIONS FOR NUCLEOTIDE DIVERSITY PIN/PIS
## Modified version - April, 2019

##  Test file -. enviroments example_sq
#load("./data/guide_pc.RData")

##### FUNCTIONS
################
##  NS SITES OF A SEQUENCE -> depends on calculate_NS_codon previously calculated.
calculate_NS_sites = function(r) {
  if(class(r) == 'DNAString') r = DNAStringSet(r)
  n_aa = width(r)/3
  s = 0
  for(n in 1:n_aa){  ##f per codon.
    codon = substr(r, 3*n-2, 3*n)
    s = s + guide$S[[codon]]
  }
  n = (3* n_aa - s)
  sites=c(s, n); names(sites) = c('S', 'N')
  return(sites)
}

## piNS in haplotypes using pairwise distances.
calculate_piNS = function(set, freqs=NULL, sites=NULL) {
  require(magrittr)
  if(length(set)<=1){
    piSN = rep(0, 3)
  }else{
    if(is.null(sites)){
      sites = sapply(set, calculate_NS_sites) %>% t ## Calculate NS sites if info is not given.
    }else{
      if(class(sites) != 'matrix') sites = as.matrix(sites)
    }
    if(is.null(freqs)) freqs = rep(1, length(set))
    piSN = 0  ## Initializing data
    for(i in 2:length(set)){
      for(j in 1:i){
        sn = c(0, 0)
        ## Retrieve s and n from each pair of codons -> codons where there is indeed a difference.
        cm = consensusMatrix(set[c(i,j)])[1:4, ]
        ix = apply(cm, 2, function(i) sum(i ==1) == 2)
        ix = which(ix > 0); ix = ceiling(ix/3) %>% unique
        for(k in ix){
          codons_k = substr(set[c(i,j)], 3*k-2, 3*k)
          sn = sn + pc[[codons_k[1]]][[codons_k[2]]] %>% unlist
        }
        ## Calculate dN and dS  in the sequences -> Average S and N between pair of sequences
        SN = apply(sites[c(i,j),], 2, mean)
        dSN = -0.75 * log(1 - 4/3 * sn/SN)
        ##  Add the value to the global sum -> Weight with absolute frequencies of  i and j seq.
        piSN = piSN +  (dSN * freqs[i] * freqs[j])
      }
    }
    n = sum(freqs)
    piSN = piSN / ((n^2 - n)/2)  ## Divide sum by number of pairs without permutations.
    pi = nt_diversity_hap(set, freqs, aln=T)  ## Calculate pi (just from sequences)
    piSN = c(pi, piSN)
  }
  names(piSN) = c('pi', 'piS', 'piN')
  return (piSN)
}

##piNS in haplotypes in differente windows-> dataframe with piNS calculations for different positions within the sequence. (uses calculate_piNS)
calculate_piNS_windows = function(set, freqs, ws=7){
  piNS = data.frame(pos=0, pi=0, piN=0, piS=0)
  ##  Actually search for those codons where there is at least 1 mutation.
  cm = consensusMatrix(set)[1:4, ]
  ix = apply(cm, 2, function(i) length(unique(i))-2)
  ix = which(ix > 0); ix = ceiling(ix/3) %>% unique
  n_aa = width(set)[1] / 3
  ##  Go in windows of size 7 aa.
  #ws = 7 ##  Make in odd
  for(j in ceiling(ws/2):(n_aa-floor(ws/2))){
    vars = substr(set, (j-floor(ws/2))*3-2, (j+floor(ws/2))*3) %>% factor
    fq = tapply(freqs, vars, sum) ## Freqs of 'local haplotype'
    codon = levels(vars) %>% DNAStringSet
    if(length(levels(vars)) == 1){
      piNS = rbind(piNS, data.frame(pos=j, pi=0, piN=0, piS=0))
    }else{
      pi = c(j, calculate_piNS(set=codon, freqs=fq) )
      piNS = rbind(piNS, pi)
    }
  }
  return(piNS[-1,])
}

###NT diversity created for haplotypes.  --> raw version  Freq = number of fragments per haplotype.
nt_diversity_hap = function(seq, freq=NULL, aln=T){ ##Input = DNAStringSet + freq in seq names or as independent vector.
  library(ape)
  if(aln==F){
    data = make_aln(seq, bin=T)
  }else{
    data = as.DNAbin(seq)
  }
  n=length(data)
  #L=max(width(seq))
  dist = as.matrix(dist.dna(data, 'raw'))
  ix = colnames(dist)
  if(is.null(freq))   freq = as.numeric(matrix(unlist(strsplit(ix, '_t')), 2, n)[2, ])
  ###pi = nucleotide diversity.
  pi = 0
  for(i in 1:length(freq)){
    pi = pi + sum(freq[i] * freq[1:i] * dist[i, 1:i])
  }
  pi = 2 * pi / (sum(freq) * (sum(freq)-1))  ##  Assuming that num of reads is larger than 2000, then n^2 ~  n(n-1)
  #return(data.frame(pi=pi, var=var))
  return(pi)
}
