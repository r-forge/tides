TidalCharacteristics <- function (	h,  		#(Water level) time series. data frame with time and h column
					h0 = h$h0, 	#Reference level, either single valued or vector with dimension corresponding to h
					h0marg = 0.3,   #[Maybe obsolete]Margin on reference level, to cope with small fluctuations in the Water level time series
					T2 = 5*60*60, 	#'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h
					filtconst = 50,	#Filtering constant for smoothing the time series
					dtMax = 15,	#maximum accepted time interval in a continuous series. Bigger time intervals are considered to be gaps
					unit = "mins",  # unit of dtMax, Tavg
					Tavg = 12.4*60) #Average period of time 
{
if (any(!is.element("time",names(h)))) stop("h must be a data frame with columns time (POSIXt) and h (numeric)")
if (class(h$time[1])[1]!="POSIXt") stop("h must be a data frame with columns time (POSIXt) and h (numeric)")
if (is.null(h0)) stop("provide reference level h0")

#Calculate high and low water levels (extrema) in water level time series
output <- extrema(h=h,h0=h0,T2 = T2,filtconst=filtconst)
HL <- output$HL
tij <- output$h

###
#determine gaps in 'continuous' data
gaps <- gapsts(tij,dtMax=dtMax)
if (!is.null(gaps)) 	{gaps$N <- tij$N[match(gaps$t1,tij$time)]	#N counts the tidal cycles
			tij$n <- findInterval(tij$time,gaps$t2)+1	#n counts the continuous series of cycles.	
} else tij$n <- 1

###
#Calculate inundation times and dry times
ITDT <- IT(tij[c("time","h")],h0=tij$h0,dtMax=dtMax)
ITs <- ITDT$IT
ITs$N <- tij$N[match(ITs$t1,tij$time)]
#When gaps are present in the time series, remove all broken cycles +- 1 cycle number to be sure, i.e. in which a gap of data exists
if (!is.null(gaps)) ITs <- subset(ITs,!is.element(N,c(gaps$N-1,gaps$N,gaps$N+1)))
  
DTs <- ITDT$DT
DTs$N <- tij$N[match(DTs$t1,tij$time)]
#When gaps are present in the time series, remove all broken cycles +- 1 cycle number to be sure
if (!is.null(gaps)) DTs <- subset(DTs,!is.element(N,gaps$N)&!is.element(N,gaps$N+1)&!is.element(N,gaps$N - 1))

###
#Calculate total inundation frequence
if (is.null(gaps)) gapstime <- 0 else gapstime <- sum(unclass(subset(gaps,t2<end)$dt))
Ncycles <- floor(unclass(difftime(min(max(tij$time,na.rm=T),end),min(tij$time,na.rm=T),unit="mins") - gapstime)/(Tavg))[1]
IF <- IF(subset(HL,HL=="H"&time<end),subset(HL,HL=="H"&time<end)$h0,N=Ncycles)


TideChars <- list(HL=HL,h=tij,gaps=gaps,IF=IF,ITs=ITs,DTs=DTs,h0 = h0,Ncycles=Ncycles)
class(TideChars) <- "Tides" 
return(TideChars)
}

print.Tides <- function(x,...){
  if (is.null(x$gaps)) cat("There were no gaps in the time series","\n") else cat("The time series consists of",max(x$gaps$n),"continuous sub-series","\n")
  cat("Time span:			", x$Ncycles, "average (tidal) cycles","\n")
  cat("Inundation frequency:		", x$IF, " (",x$IF*x$Ncycles/100,"inundations during time span)","\n")
  cat("Average inundation height:	", mean(x$HL$h-x$HL$h0),"\n")
  cat("Average inundation time:	", mean(x$ITs$dt),x$unit,"\n")
  cat("Average dry time:		", mean(x$DTs$dt),x$unit,"\n")

}

plot.Tides <- function(x,...){
  plot(x$h$time,x$h$h,type="l",ylab="waterlevel")
  lines(x$h$time,x$h$h0)
  points(x$HL$time,x$HL$h,col="red",pch=20)
  points(x$HL$time[x$HL$HL=="L"],x$HL$h[x$HL$HL=="L"],col="blue",pch=20)
}


extrema <- function(h,            #(Water level) time series. data frame with time and h column
                    h0,           #Reference level, either single valued or vector with dimension corresponding to h
                    h0marg = 0.3, #Margin on reference level, to cope with small fluctuations in the (water level) time series
                    T2 = 5*60*60, #'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h
                    filtconst = 50 #Filtering constant for smoothing the time series
)
{

h$h0 <- h0

#set all levels < h0 equal to h0
#take a margin (default is 0.3 cm), to cope with small fluctuations of water level
h$h[floor(h$h)<=h$h0+h0marg] <- h$h0[floor(h$h)<=h$h0+h0marg]

#filter tij data, running average of filtconst (default = 3) succesive datapoints,to
#remove small fluctuations
h$ho <- h$h
h$hfilt <- filter(h$h,rep(1/filtconst,filtconst))

#remove missing values due to filtering
h <- h[!is.na(h$hfilt),]

#Add missing values due to filtering
#h$h[1:floor(filtconst/2)] <- h$ho[1:floor(filtconst/2)]
#h$h[(length(h$h)-floor(filtconst/2)+1):length(h$h)] <- h$ho[(length(h$h)-floor(filtconst/2)+1):length(h$h)]

#Useful matrix to swith between [H,L] or [T,F] representation of high and low phase of ts
HLTF <- data.frame(HL = c("H","L"), TF = c(TRUE,FALSE))

#Here the core thing happens
#If h[t+T2] < h[t] & h[t-T2] < h[t] then high tide, else low tide
h$TF <- (h$hfilt>approx(x=h$time,y=h$hfilt,xout= pmin(h$time+T2,h$time[length(h$time)]))$y)+h0marg&(h$hfilt>approx(x=h$time,y=h$hfilt,xout=pmax(h$time[1],h$time-T2))$y+h0marg)
h$HL <- HLTF$HL[match(h$TF,HLTF$TF)]

#Give every high and low tide phase a number
h$N <- 0
h$N[2:(length(h$time))] <- 1*(h$HL[1:(length(h$time)-1)] != h$HL[2:(length(h$time))])
h$N[1] <- 1
h$N <- cumsum(h$N)

#find all maxima within each high and low phase
max <- as.data.frame(tapply(h$h,h$N,max))
names(max) <- "max"
min <- as.data.frame(tapply(h$h,h$N,min))
names(min) <- "min"

h$max <- max$max[h$N]
h$min <- min$min[h$N]

#Tom check this
#This way, sometimes two extrema (with the same water level as the extremum) are counted
h$HLval <- 0
h$HLval <- (h$max==h$h & h$HL=="H")*h$h + (h$min==h$h & h$HL=="L")*h$h

#select only High and Low waters
HL <- h[h$HLval != 0,]

#select only different extrema, take first
HL$eq <- c(F,HL$HLval[1:(length(HL$time)-1)] == HL$HLval[2:length(HL$time)])
HL <- HL[!HL$eq,]
return(list(HL = HL[c("time","h","HL","h0")], #Data frame with extrema
              h = h[c("time","h","h0","HL","N")]) #original water level data frame with additional attributes
              )
}

gapsts <- function(ts,            # data-frame with time column
                    dtMax,        # maximum accepted time interval in a continuous series
                    unit = "mins"  # unit of dtMax
                    )
{
#Select gaps > dtMax in a timeseries ts
timediffs <- difftime(ts$time[1:(length(ts$time)-1)], ts$time[2:(length(ts$time))],units=unit)
if (!any(timediffs < - dtMax)) return(NULL)
gaps <- ts$time[c(timediffs < -dtMax,F)]
gaps <- data.frame(t1 = gaps)
gaps$t2 <- ts$time[match(gaps$t1,ts$time)  + 1]
gaps$n <- 1:dim(gaps)[1]
gaps$dt <- difftime(gaps$t2,gaps$t1,units=unit)
return(gaps) #Data frame with the initial time, end time and time difference (unit = unit) of each interval > dtMax
}

IT <- function(h,             #Water level time series. data frame with time and h column
                h0,           #Reference level, either single valued or vector with same length as h
                h0marg = 0.3, #Margin on reference level, to cope with small fluctuations in the Water level time series
                dtMax = 15,   #Maximum time interval in continuous water level series
                unit = "mins"  #Unit of dtMax and output time intervals
                )
{

dry <- subset(h[c("time","h")],h<=(h0 + h0marg))
if (dim(dry)[1] == 0) {		
#If the site never falls dry, inundation time equals the time of the time series
  IT <- data.frame(t1 = h$time[1],t2 = h$time[length(h$time)], dt = difftime(h$time[length(h$time)],h$time[1],units = unit))
  DT <- NULL
} else 
{
wet <- subset(h[c("time","h")],h>(h0 + h0marg))	#dry time = 'inverse' of inundation time
if (dim(wet)[1] == 0) {	
#If the site is never inundated, dry time equals the time of the time series
  IT <- NULL	
  DT <- data.frame(t1 = h$time[1],t2 = h$time[length(h$time)], dt = difftime(h$time[length(h$time)],h$time[1],units = unit))	
} else
{
if (wet$time[1] > h$time[1]) {
			wet <- rbind(h[1,],wet)
			wet$time[1] <- wet$time[1] - unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (wet$time[length(wet$time)] < h$time[length(h$time)]){
			wet <- rbind(wet,h[length(h$time),])
			wet$time[length(wet$time)] <- wet$time[length(wet$time)] + unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (dry$time[1] > h$time[1]) {
			dry <- rbind(h[1,],dry)
			dry$time[1] <- dry$time[1] - unclass(difftime(h$time[2],h$time[1],units="secs"))
			}
if (dry$time[length(dry$time)] < h$time[length(h$time)]) {
			dry <- rbind(dry,h[length(h$time),])
			dry$time[length(dry$time)] <- dry$time[length(dry$time)] + unclass(difftime(h$time[2],h$time[1],units="secs"))
			}

IT <- gapsts(dry,dtMax,unit=unit)
DT <- gapsts(wet,dtMax,unit=unit)

}
}

return(list(IT = IT,     #Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of inundation
            DT = DT))    #Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of dry time
}

IF <- function(H,                  #High water levels. Data Frame with column h
                h0,                #Reference level for which IF has to be calculated, either single valued or array of length = length(H[1,])
                N  = length(H[,1]) #number of cycles in time series, equals the number of high water levels when these are complete (= default value)
                )
{
IF <- length(subset(H[c("h")],h>h0)[,1])/N
return(IF*100)               #Inundation frequence [%] at (varying) reference level h0)
}
