#######
## Miscelaneous functions.
#######

#########
## Finding frame
########

##  Determine frame and get in frame.
#####  Function to find correct reading frame
determine_frame = function(sq){
  require('Biostrings')
  if(length(sq) > 1) sq = sq[1]
  ln = width(sq)
  ##  First:  Try all 6 possible options.
  fw = suppressWarnings( sapply(1:3, function(i) translate(subseq(sq, i, ln)))  )
  rv = suppressWarnings( sapply(1:3, function(i) translate(subseq(reverseComplement(sq), i, ln))) )
  stop_codons = sapply(c(fw, rv), letterFrequency, letters=c('*'))
  frame_ix = which(stop_codons == 0)
  if(length(frame_ix) != 1) stop('None or various possible reading frames')
  ##  Choose correct sense.
  sense = function(sq) return(sq)
  if(frame_ix > 3) sense = function(sq) return(Biostrings::reverseComplement(sq))
  assign('get_inframe.sense', sense, envir = .GlobalEnv)

  ##  Choose the correct start and end point for the sequence.
  x = frame_ix %% 3
  y = ln - (4 - which( (ln-2):ln %% 3 == x ) )
  if(x == 0) x = 3  ##Starting position.
  assign('get_inframe.xy', c(start = x,stop = y), envir = .GlobalEnv)
  return(get_inframe.xy)
}

##  New get_inframe function.

get_inframe = function(sq){
  sq2 = get_inframe.sense(sq)
  xy = get_inframe.xy
  return( subseq(sq2, xy[1], xy[2]) )
}



#######
## Others
#######
