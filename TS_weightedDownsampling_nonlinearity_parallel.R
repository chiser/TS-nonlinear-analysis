################## An R-script to import strokelitude data and to plot it in several ways

## if you did not do so yet, set the right working directory
setwd("C:/Users/LocalAdmin/Desktop")

## source the script with the functions needed for analysis
source("StrokePrepFunctions.R")
require("nonlinearTseries")
require("beepr")
#filename = id_table1[i]
Prediction<- NULL
Prediction2<- NULL
Params<-NULL
Observed<-NULL
Observed2<-NULL
# esto se va a incluir en futuras versiones del paquete, es util para utilizar
# vectores en el argumento prediction.step... recomendado!
nonLinearPredict = Vectorize(nonLinearPrediction, vectorize.args = "prediction.step")
nonLinearPredict2 = Vectorize(nonLinearPredict, vectorize.args = "radius")

message("Enter the prefix for the saved files")
Textname<- scan(n=1,what= character())

## do the sampling by taking bins of time (instead of samples). By Christian Rohrsen
#print("Please enter the time between samples in seconds")
#samplePeriod <- 0.05   #scan(n=1)   # set the time of sampling in seconds

print("Please enter the number of flies")
Nflies <- scan(n=1) 



knnfunction<- function(Prediction,Prediction2,Observed,Observed2,nonLinearPredict2) {

##### read the data with the corresponding function
beep(2)
flyTraces <- flyDataImport()

plot(diff(flyTraces$Time[1:5000]))

print("Please enter the cutting starting point")
StartingAnalysis <- scan(n=1) 

flyTraces<-flyTraces[StartingAnalysis:length(flyTraces$Time),1:3]

plot(diff(flyTraces$Time))

print("Do you want to delete more pieces (y/n)?")
Answer <- scan(what="character",n=1) 

if(Answer=="y"){
  print("Start")
  ChangeStart <- scan(n=1)
  
  print("End")
  ChangeEnd <- scan(n=1)
  
  flyTraces<-flyTraces[ChangeStart:ChangeEnd,1:3]
}

plot(diff(flyTraces$Time))

flyTracesFiltered <- flyDataFilter(flyTraces, frequency = .12, order = 3)

flyTracesFiltered$group_num <- floor(flyTracesFiltered$Time/0.05)

flyTracesFiltered2<-NULL


### Downsampling but weighting according to the measured time within the bin

flyTracesFiltered$weight <- (abs(flyTracesFiltered$Time-(flyTracesFiltered$group_num*0.05+0.025)))^(-1)

flyTracesFiltered$Norm<-ave(flyTracesFiltered$weight,flyTracesFiltered$group_num,FUN=function(x) x/sum(x))

flyTracesFiltered$Time2<- flyTracesFiltered$Time*flyTracesFiltered$Norm

flyTracesFiltered$Right2<- flyTracesFiltered$Right*flyTracesFiltered$Norm

flyTracesFiltered$Left2<- flyTracesFiltered$Left*flyTracesFiltered$Norm

flyTracesFiltered2$Time<- tapply(flyTracesFiltered$Time2, flyTracesFiltered$group_num, sum)

flyTracesFiltered2$Right<- tapply(flyTracesFiltered$Right2, flyTracesFiltered$group_num, sum)

flyTracesFiltered2$Left<- tapply(flyTracesFiltered$Left2, flyTracesFiltered$group_num, sum)

flyTracesFiltered2$Trace = flyTracesFiltered2$Right-flyTracesFiltered2$Left

use.ts= flyTracesFiltered2$Trace


### Downsampling without any weighting

#TracesFiltered2<-NULL

#TracesFiltered2$Right = aggregate(flyTracesFiltered$Right,
#                               list(group_num=flyTracesFiltered$group_num), mean)
#TracesFiltered2$Left = aggregate(flyTracesFiltered$Left,
#                                     list(group_num=flyTracesFiltered$group_num), mean)
#TracesFiltered2$Time = aggregate(flyTracesFiltered$Time,
#                                     list(group_num=flyTracesFiltered$group_num), mean)

#TracesFiltered2$Trace = TracesFiltered2$Right$x-TracesFiltered2$Left$x



LAG <- timeLag(use.ts, technique = "ami",selection.method = "first.e.decay",lag.max = 20, do.plot = TRUE)


# selecting the embedding dim ---------------------------------------------

emb.dim = estimateEmbeddingDim(use.ts, threshold = 0.95,
                               time.lag = LAG, max.embedding.dim = 25)


#spaceTime = spaceTimePlot(time.series = use.ts,time.lag = LAG,
#                          embedding.dim = emb.dim,
#                          number.time.steps = 40,type="b")

#message("please enter the theiler window according to the graph")
#theilerw = scan(n=1,what= numeric())
#theilerw = 30

#corr.dim =  corrDim(use.ts, min.embedding.dim = emb.dim, 
#                    max.embedding.dim = emb.dim + 7,
#                    time.lag = LAG,
#                    min.radius = 0.04, max.radius = (max(use.ts)-min(use.ts))/2, 
#                    n.points.radius = 45, 
#                    theiler.window = theilerw,
#                    do.plot = TRUE)
#beep(1)
#message("please enter the parameters for the estimation of the correlation dimension")## TODO: take a corr.dim$corr.matrix and find the max slope and so it is automatically
#regression.range.min = scan(n=1,what= numeric())
#message("please enter the parameters for the estimation of the correlation dimension")
#regression.range.max = scan(n=1,what= numeric())
#message("please enter the parameters for the estimation of the correlation dimension")
#use.embeddings.min = corr.dim$embedding.dims[5]
#message("please enter the parameters for the estimation of the correlation dimension")
#use.embeddings.max = corr.dim$embedding.dims[8]


#D<- estimate(corr.dim,regression.range = c(regression.range.min,regression.range.max),use.embeddings = use.embeddings.min:use.embeddings.max,do.plot = T)



#detrend = dfa(use.ts, 
#              window.size.range = c(10, 400),
#              npoints = 20,do.plot=TRUE)     
#detr.fluct<-estimate(detrend,regression.range = c(30,300))

#ml<-maxLyapunov(use.ts, min.embedding.dim = 2,                                        #from nonlinearTseries
#                max.embedding.dim = emb.dim+4, time.lag = LAG, radius=1,
#                theiler.window = theilerw, min.neighs = 5, min.ref.points = 500,
#                max.time.steps = 10,
#                do.plot = TRUE)

#plot(ml)
#ml.estimation = estimate(ml,regression.range = c(0,15),
#                         use.embeddings=2:emb.dim+4,
#                         do.plot = TRUE)

#surrogate data test

#st = surrogateTest(use.ts,significance = 0.05,one.sided = F,
#                   FUN = timeAsymmetry, do.plot=T)

#sample entropy

#se = sampleEntropy(corr.dim, do.plot = T)
#se.est = estimate(se, do.plot = T,
#                  regression.range = c(0,20))
#EstimatedEntropy<- mean(se.est)

save.image()

for (oo in 1:8){
  beep(1)
  message("please enter the radius")
  radius = scan(n=1,what= numeric())
  
  # obtenemos predicciones para las proximas 500 (prediction.step) muestras basandonos
  #en las 4000 primeras muestras de x
  prediction = nonLinearPredict2 (use.ts[1:round(length(use.ts)/2)], embedding.dim = emb.dim, time.lag = LAG, radius = c(0.01,radius),
                                  radius.increment = 0.02, prediction.step = 1:40)
  # comparamos... Fijate como la prediccion va empeorando a medida que pasa el tiempo. Esto
  # es tipico de sistemas no lineales
  
  plot(use.ts[(round(length(use.ts)/2)+1):(round(length(use.ts)/2)+40)], type = "l", ylim = range(use.ts[(round(length(use.ts)/2)+1):(round(length(use.ts)/2)+40)], Prediction))
  lines(prediction[,1], type = "l", col = 2)
  lines(prediction[,2], type = "l", col = 3)

  
  beep(2)
  message("Is the prediction fine?y/n")
  PredictionAnswer = scan(n=1,what= character())
  
  if(PredictionAnswer== "y") break
  
}

for (oo in 1:8){
  beep(1)
  message("please enter the radius")
  radius2 = scan(n=1,what= numeric())
  
  # obtenemos predicciones para las proximas 500 (prediction.step) muestras basandonos
  #en las 4000 primeras muestras de x
  prediction2 = nonLinearPredict2 (use.ts[1:(length(use.ts)-2000)], embedding.dim = emb.dim, time.lag = LAG, radius = c(0.01,radius2),
                                   radius.increment = 0.02, prediction.step = 1:40)
  # comparamos... Fijate como la prediccion va empeorando a medida que pasa el tiempo. Esto
  # es tipico de sistemas no lineales
  
  plot(use.ts[((length(use.ts)-2000)+1):((length(use.ts)-2000)+40)], type = "l", ylim = range(use.ts[((length(use.ts)-2000)+1):((length(use.ts)-2000)+40)], Prediction2))
  lines(prediction[,1], type = "l", col = 2)
  lines(prediction[,2], type = "l", col = 3)

  
  beep(2)
  message("Is the prediction fine?y/n")
  PredictionAnswer = scan(n=1,what= character())
  
  if(PredictionAnswer== "y") break
  
}



#manytests<-nonlinearityTest (time.series=use.ts, verbose = TRUE)


#denoised<-nonLinearNoiseReduction(time.series=use.ts, embedding.dim=emb.dim, radius=0.1)
#plot(denoised,type="l")

#print(D,rqa.res,mean(se.est),ml.estimation)   


#dev.off()

Observed<-c(Observed,use.ts[(round(length(use.ts)/2)+1):(round(length(use.ts)/2)+40)])
Observed2<-c(Observed2,use.ts[((length(use.ts)-2000)+1):((length(use.ts)-2000)+40)])
Prediction<- cbind(Prediction,prediction)
Prediction2<- cbind(Prediction2,prediction2)

knn<-list(Observed,Observed2,Prediction,Prediction2)

return(knn)

}


see<- mclapply.hack(Nflies,knnfunction(Prediction,Prediction2,Observed,Observed2,nonLinearPredict2))

write.table(Params, paste(Textname,"Spontaneity.txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE)
write.table(Prediction, paste(Textname,"Pred.txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE)
write.table(Prediction2, paste(Textname,"Pred2.txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE)
write.table(Tests, paste(Textname,"Tests.txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE)
write.table(Observed, paste(Textname,"Observed.txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE)
write.table(Observed2, paste(Textname,"Observed2.txt",sep=""), sep="\t", row.names = FALSE,col.names = FALSE)

beep(3)