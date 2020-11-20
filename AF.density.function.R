# AF density function
# Author Efstratios I. Charitos, MD, PhD efstratios.charitos@gmail.com
# Date: 20-Nov-2020 v.4.2
# Code optimization and speed improvements Joao Moteiro, PhD joao.v.monteiro@medtronic.com

#data is a vector of daily AF minutes
#timeunits = how many time units (data) are in one day : If data is in days then timeunits=24; if data is in minutes then timeunits=1440
#limit=T  hardcodes a density of 1 in cases of >1 results (high burdens)
#the switch expand.vector.resolution.minute.level=F expands the input vector of daily AF minutes to 1 minute intervals. This handles single day episodes better (for an imput vector of daily AF minutes)



AF.density.v0.4.2 <- function(data, minim=0.001, timeunits=1440, limit=T, expand.vector.resolution.minute.level=T){
  
  library(RcppRoll)
  
  
  f2 <- function(v, minim, maxim, b, timeunits) {
    tx <- seq(minim, maxim, b)
    min_n_units <- ceiling(tx*sum(v)) #minimum number of 1s needed to reach desired prop
    maxM <- length(min_n_units)
    k <- rep(Inf, maxM)
    if(timeunits == 1) {
      
      rle_v <- rle(v)
      cont_blocks <- rle_v$lengths[rle_v$values == 1]   #length of blocks of contiguous 1s
      max_n_blocks <- length(cont_blocks)
      #find the lengths fit within a block of contiguous ones 
      max_cont_blocks <- max(cont_blocks)
      k[min_n_units <= max_cont_blocks] <- min_n_units[min_n_units <= max_cont_blocks]   
      
      zero_spaces <- rle_v$lengths[rle_v$values == 0]
      if(rle_v$values[1] == 0) zero_spaces <- zero_spaces[-1]   
      if(rev(rle_v$values)[1] == 0) zero_spaces <- zero_spaces[-length(zero_spaces)]
      
      if(!is.finite(k[maxM])) {
        nblocks <- 2
        
        while(nblocks <= max_n_blocks) {
          len <- roll_sum(cont_blocks, nblocks)
          pad <- roll_sum(zero_spaces, nblocks - 1)
          
          pos_to_look <- which(k > max_cont_blocks & min_n_units <= max(len))
          
          for(i in pos_to_look) {
            k[i] <- min(k[i], min_n_units[i] + min(pad[len-min_n_units[i]>= 0]))
          }
          nblocks <- nblocks + 1
        }
      }
    } else {
      length_int <- 1
      while (!is.finite(k[maxM])) {
        k[max(roll_sum(v,length_int)) >= min_n_units & !is.finite(k)] <- length_int
        length_int <- length_int + 1
      }
    }
    
    return(k/length(v))
  } 
  
  
  if ( expand.vector.resolution.minute.level==T) {
    timeunits=1
    #Vector expansion to minute level
    data <- rep(rep(c(1, 0), length(data)), c(rbind(data, 1440-data)))
  }
  
  
  if (sum(data)==0) return(0) else { #CHECK
    
    maxim=1-minim
    b=minim
    
    data.test<-data.frame(day=1:length(data), afmin=data)
    burden<-sum(data.test$afmin)/timeunits/length(data.test$afmin)
    data.test$daily.burd<-data.test$afmin#/sum(data.test$afmin)
    
    raw.den<-(2*(sum(abs(seq(minim, maxim, b)-(f2(data.test$daily.burd, minim=minim, maxim=maxim, b=b, timeunits=timeunits)))*b))/(1-burden))
    
    
    density <- raw.den
    if (limit==T & raw.den>1) density=1
    #if (is.na(raw.den)==TRUE) density=0
    
    return(density)
  }
}
