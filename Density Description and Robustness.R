# AF density Description and Robustness
# Author Efstratios I. Charitos, MD, PhD efstratios.charitos@gmail.com
# Date: 22-Nov-2019 v.2


####
# Example Patients

pt.A<-c(rep(1440, 80), rep(0,220))
pt.B<-c(rep(0,110),rep(1440, 80), rep(0,110))

pt.C<-c(rep(0,59), seq(0,1440, 9), rep(0,80))
pt.D<-(rev(pt.C))

pt.E<-c(rep(0,109), runif(81,min = 0, max = 1440), rep(0,110))
pt.F<-c(rev(pt.E))

pt.Z<-c(0, 0, 0, 0, 360, 0, 0, 0, 0)
AF.density(pt.Z)

pt.Z.minute.level.resolution<-c(rep(0,4*1440), rep(1,360), rep(0,4*1440))
AF.density(pt.Z.minute.level.resolution, timeunits = 1)


pt.Z.minute.level.resolution.7d<-c(rep(0,3*1440), rep(1,5), rep(0,3*1440))
pt.Z.minute.level.resolution.14d<-c(rep(0,6*1440), rep(1,5), rep(0,7*1440))
pt.Z.minute.level.resolution.30d<-c(rep(0,14*1440), rep(1,5), rep(0,15*1440))

AF.density(pt.Z.minute.level.resolution.7d, timeunits = 1)
AF.density(pt.Z.minute.level.resolution.14d, timeunits = 1)
AF.density(pt.Z.minute.level.resolution.30d, timeunits = 1)



#Function


plot.compass.index<-function(x,roll.mean.n=2,minim=10,maxim=90,b=20,adjust=40){

	
#F2 index subfunction
f2.index<-function(minim=0.01, maxim=0.99, b=0.01, data){
			f.my<-function(x, max) {
				
			s <- 0
			len <- Inf
			start <- 1
			j <- 1
			for (i in seq_along(x)) {
				s <- s + x[i]
				while (s >= max) {
					if (i-j+1 < len) {
						len <- i-j+1
						start <- j}
					s <- s - x[j]
					j <- j + 1}}
					
		list(s=start, l=len)}
		start<-NULL;len<-NULL;end<-NULL
		for (m in seq(from=minim, to=maxim, by=b)){
		start<-c(start,f.my(data,m)$s)
		len<-c(len,f.my(data,m)$l)	} 
		final<-data.frame(start=start, len=len, end=start+len)
		return(final)}
#F2 function
f2<-function(minim=0.01, maxim=0.99, b=0.01,data){
	#data=burden per day
	f<-function(x, max) {
		
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
				j <- j + 1}}
		list(start=start, length=len)
		len}
		k<-NULL
	for (m in seq(from=minim, to=maxim, by=b)){ k<-c(k,f(data,m))} 
	k.per<-k/length(data)
	return(k.per)}

#AF Density
AF.density<-function(data, minim=0.01, timeunits=1440){
	
	if (sum(data)==0) return(0) else { #CHECK
		
		#data is a vector of daily AF minutes
		
		maxim=1-minim
		b=minim
		
		data.test<-data.frame(day=1:length(data), afmin=data)
		burden<-sum(data.test$afmin)/timeunits/length(data.test$afmin)
		data.test$daily.burd<-data.test$afmin/sum(data.test$afmin)
		
		raw.den<-(2*(sum(abs(seq(minim, maxim, b)-(f2(data.test$daily.burd, minim=minim, maxim=maxim, b=b)))*b))/(1-burden))
		
		density<-raw.den
		if (raw.den>1) density=1
		if (is.na(raw.den)==TRUE) density=0
		
		return(density)
	}
}

	
	burden.pat<-sum(x)/(1440*length(x))
	tot.mon.days.pat<-length(x)
	frame<-f2.index(data=(x/sum(x)))
	
	
par(mfrow=c(1,3))
par(mar=c(5.1,5.1,4.1,2.1))
	plot(1:length(x), x, xlim=c(1,length(x)),type="h", ylim=c(5,1441), xlab="Monitored Days", ylab="Daily min in AF",col="gray1", cex.lab=1.25, cex.axis=1.2, main="Rhythm History",cex.main=1.8, bty="n")
	legend("topright", bty="n",c(paste("AF density:",round(AF.density(x),2))  ,paste("AF burden:",round(burden.pat,2))),cex=1)
	#Starting index
	for (i in seq(minim,maxim,b))  {
		lines(c(frame$start[i],frame$end[i]), c(i*1440/100, i*1440/100), col="red",lwd=7)
		text(x=frame$start[i]+adjust,y=(i*1440/100)+50, paste(i,"%: ",frame$len[i]," days", sep="" ),cex=1.2, col="red")}


plot(f2(data=x/sum(x)),seq(0.01,0.99,0.01),  lwd=5,type="l", xlim=c(0,1), ylim=c(0,1), ylab="Proportion of burden", xlab="Proportion of time", lty=3, col="darkblue", cex.lab=1.25,cex.axis=1.2, bty="n", main="Observed",cex.main=1.8)
curve(2*x/2, from=0,to=1,add=T, lwd=4)
polygon(y=c(0,seq(0.01,0.99,0.01),1,0), x=c(0, (f2(data=x/sum(x))),1,0), col=rgb(149,149,219,255,maxColorValue=255), border=NA)
lines((f2(data=x/sum(x))),seq(0.01,0.99,0.01),  lwd=5,type="l", xlim=c(0,1), ylim=c(0,1), ylab="Proportion of burden", xlab="Proportion of time", lty=3, col="darkblue", cex.lab=1.25,cex.axis=1.2)
curve(2*x/2, from=0,to=1,add=T, lwd=4)
points(f2(data=x/sum(x))[c(10,30,50,70,90)], c(.1,.3,.5,.7,.9), pch=20, cex=2)
blue<-sum(0.01*abs(seq(0.01,0.99,0.01)-f2(data=x/sum(x))))
legend("bottomright", paste("Blue area: ", round(blue,2)), bty="n", cex=1.5)

#Add the max den lines
#max.den<-c(rep(1440, burden.pat*length(x)), rep(0,(1-burden.pat)*length(x)))
#max.den<-((max.den)/sum(max.den))
#plot(f2(max.den,minim=0.01, maxim=0.99,b=0.01), seq(0.01,0.99,0.01), lty=1, col="darkgreen", lwd=5,type="l", xlim=c(0,1), ylim=c(0,1), ylab="Proportion of burden", xlab="Proportion of time",   cex.lab=1.25,cex.axis=1.2, bty="n", main="Maximum Density",cex.main=1.8)

plot(0,0, lty=1, col="darkgreen", lwd=5,type="l", xlim=c(0,1), ylim=c(0,1), ylab="Proportion of burden", xlab="Proportion of time",   cex.lab=1.25,cex.axis=1.2, bty="n", main="Maximum Density",cex.main=1.8)
lines(c(0,burden.pat), c(0,1), lty=3, lwd=5, col="darkgreen")
curve(2*x/2, from=0,to=1,add=T, lwd=4)

polygon(y=c(0,1,1), x=c(0,1,burden.pat), col=rgb(120,190,120,255,maxColorValue=255), border=NA)
green<-(1-burden.pat)/2
legend("bottomright", paste("Green area: ", round(green,2)), bty="n", cex=1.5)
}










#Testing

plot.compass.index(pt.A)
plot.compass.index(pt.B)
plot.compass.index(pt.C)
plot.compass.index(pt.D)
plot.compass.index(pt.E)
plot.compass.index(pt.F)
