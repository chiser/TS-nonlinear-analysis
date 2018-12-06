################## An R-script to import strokelitude data and to plot it in several ways

## if you did not do so yet, set the right working directory
setwd("C:/Users/LocalAdmin/Desktop/scripts/repositories/StrokelitudeScripts")

## source the script with the functions needed for analysis
source("StrokePrepFunctions.R")
require("nonlinearTseries")

setwd("C:/Users/LocalAdmin/Desktop/data/Nonlinearity results/c105c232_strokelitude")

  ##### read the data with the corresponding function
  flyTraces <- flyDataImport()
  
  flyTraces$Trace<- flyTraces$Right-flyTraces$Left
  
  plot(flyTraces$Trace,type="l")
  
  plot(diff(flyTraces$Time))

  flyTracesFiltered <- flyDataFilter(flyTraces, frequency = .12, order = 3)
  
  plot(flyTracesFiltered$Trace,type="l")
  
  LAG <- timeLag(flyTracesFiltered$Trace, technique = "ami",selection.method = "first.e.decay",lag.max = 20, do.plot = TRUE)  
  
  use.ts<-diff(flyTracesFiltered$Trace)
  
  LAG <- timeLag(use.ts, technique = "ami",selection.method = "first.e.decay",lag.max = 20, do.plot = TRUE)
  
  
  # selecting the embedding dim ---------------------------------------------
  
  emb.dim = estimateEmbeddingDim(use.ts, threshold = 0.95,
                                 time.lag = LAG, max.embedding.dim = 25)
  
  
  spaceTime = spaceTimePlot(time.series = use.ts,time.lag = LAG,
                            embedding.dim = emb.dim,
                            number.time.steps = 40,type="b")
  

  theilerw = 30
  
  corr.dim =  corrDim(use.ts, min.embedding.dim = emb.dim, 
                      max.embedding.dim = emb.dim + 7,
                      time.lag = LAG,
                      min.radius = 0.04, max.radius = (max(use.ts)-min(use.ts))/2, 
                      n.points.radius = 45, 
                      theiler.window = theilerw,
                      do.plot = TRUE)

  
  
  D<- estimate(corr.dim,regression.range = c(regression.range.min,regression.range.max),use.embeddings = use.embeddings.min:use.embeddings.max,do.plot = T)
  

  
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
    

