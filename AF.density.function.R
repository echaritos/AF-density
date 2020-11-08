# AF density function
# Author Efstratios I. Charitos, MD, PhD efstratios.charitos@gmail.com
# Date: 22-Nov-2019 v.4 BETA

#data is a vector of daily AF minutes
#timeunits = how many time units (data) are in one day : If data is in days then timeunits=24; if data is in minutes then timeunits=1440
#limit=T  hardcodes a density of 1 in cases of >1 results (high burdens)
#the switch expand.vector.resolution.minute.level=F expands the input vector of daily AF minutes to 1 minute intervals.

AF.density.temp<-function(data, minim=0.001, timeunits=1440, limit=T, expand.vector.resolution.minute.level=F){

	if ( expand.vector.resolution.minute.level==T) {
	timeunits=1
	#Vector expansion to minute level
	data <- rep(rep(c(1, 0), length(data)), c(rbind(data, 1440-data)))
}
	

	
	
	#f2 function
	f2<-function(minim=0.001, maxim=0.999, b=0.001,data){
		
		
		f<-function(x, max) {
			
			#rolly  
			s <- 0
			len <- Inf
			start <- 1
			j <- 1
			for (i in seq_along(x)) {
				s <- s + x[i]
				while (s >= max) {
					if (i-j+1 < len) {
						len <- i-j+1
						start <- j
					}
					s <- s - x[j]
					j <- j + 1
				}
			}
			
			list(start=start, length=len)
			
			len
		}
		
		k<-NULL
		for (m in seq(from=minim, to=maxim, by=b)){ k<-c(k,f(data,m))} 
		k.per<-k/length(data)
		return(k.per)
		
	}
	
	if (sum(data)==0) return(0) else { #CHECK
	
	maxim=1-minim
	b=minim
	
	data.test<-data.frame(day=1:length(data), afmin=data)
	burden<-sum(data.test$afmin)/timeunits/length(data.test$afmin)
	data.test$daily.burd<-data.test$afmin/sum(data.test$afmin)
	
	raw.den<-(2*(sum(abs(seq(minim, maxim, b)-(f2(data.test$daily.burd, minim=minim, maxim=maxim, b=b)))*b))/(1-burden))
  
  density<-raw.den
  if (limit==T & raw.den>1) density=1
	#if (is.na(raw.den)==TRUE) density=0

  return(density)
}
}





