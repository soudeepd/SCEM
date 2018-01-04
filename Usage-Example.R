#-------------------------------------------------------------
# Preamble
#-------------------------------------------------------------

# Set the working directory 

# Install the github package
# devtools::install_github("soudeepd/SCEM") # not yet exported as a package

#-------------------------------------------------------------
# Read all the data files
#-------------------------------------------------------------
oxy <- read.csv("OxygenData3.csv")
names(oxy)[3] = "ID"
sheep <- subset(oxy,species=="Ovis")
armenia = split(sheep,f = sheep$ID)

balasse <- read.csv("BalasseData.csv")
colnames(balasse) <- c("ID","distance","oxygen")
rousay = split(balasse,f = balasse$ID)

newdata <- read.csv("newdata-balasse.csv")
bordusani = split(newdata,f = newdata$ID)


#-------------------------------------------------------------
# Implement the SCEM and Cosine method on the Armenia data
#-------------------------------------------------------------
results = SCEM(armenia,optimumbw = -0.33)
cosine = makeFits(armenia)  
rownames(cosine) = 1:nrow(cosine)

results$results$Cosine = cosine$birth
dif = abs(results$results$Cosine - results$results$Season)
results$results$Difference = apply(cbind(dif,(1-dif)),1,min)
results$results$R = cosine$Pearson
rrr.a = results$results[order(results$results$Cluster),]

# Storing the birth seasonality estimates from two methods
c.armenia = results$results$Cosine
s.armenia = results$results$Season

# Plotting the graphs of the clusters for Armenia data
gnum = max(results$results$Cluster)
dev.off()
par(mar=c(1.5,1.5,1,1), mfrow=c(2,2),oma = c(4, 4, 1, 1))
for(i in 1:gnum){
  if(length(results$groups[[i]])>1){
    cc = as.character()
    plot(c(0,40),c(-14,0),type = "n",axes = T,xlab = "",ylab = "",xlim = rev(range(c(0,40))))
    for(j in 1:length(results$groups[[i]])){
      tt = armenia[[results$groups[[i]][j]]]
      lines(tt$distance,tt$oxygen,lty = j)
      cc = c(cc,paste(tt$lab[1],tt$ID[1],sep = ""))
    }
    legend("topleft",legend = cc,lty = c(1:j),lwd = 2)
  }
}
title(xlab = "distance from CEJ (mm)",ylab = expression(delta^18*O[VPDB]*' (per mil)'),outer = TRUE, line = 2)


#-------------------------------------------------------------
# Implement the SCEM and Cosine method on the Bordusani data
#-------------------------------------------------------------
results = SCEM(bordusani,optimumbw = -0.33)
cosine = makeFits(bordusani)  
rownames(cosine) = 1:nrow(cosine)

results$results$Cosine = cosine$birth
dif = abs(results$results$Cosine - results$results$Season)
results$results$Difference = apply(cbind(dif,(1-dif)),1,min)
results$results$R = cosine$Pearson
rrr.b = results$results[order(results$results$Cluster),]

# Storing the birth seasonality estimates from two methods
c.bordusani = results$results$Cosine
s.bordusani = results$results$Season

# Plotting the graphs of the clusters for Armenia data
gnum = max(results$results$Cluster)
dev.off()
par(mar=c(1.5,1.5,1,1), mfrow=c(1,3),oma = c(4, 4, 1, 1))
for(i in 1:gnum){
  if(length(results$groups[[i]])>1){
    cc = as.character()
    plot(c(0,40),c(-14,0),type = "n",axes = T,xlab = "",ylab = "",xlim = rev(range(c(0,40))))
    for(j in 1:length(results$groups[[i]])){
      tt = bordusani[[results$groups[[i]][j]]]
      lines(tt$distance,tt$oxygen,lty = j)
      cc = c(cc,paste(tt$lab[1],tt$ID[1],sep = ""))
    }
    legend("topleft",legend = cc,lty = c(1:j),lwd = 2)
  }
}
title(xlab = "distance from CEJ (mm)",ylab = expression(delta^18*O[VPDB]*' (per mil)'),outer = TRUE, line = 2)


#-------------------------------------------------------------
# Implement the SCEM and Cosine method on the Rousay data
#-------------------------------------------------------------
results = SCEM(rousay,optimumbw = -0.33)
cosine = makeFits(rousay)  
rownames(cosine) = 1:nrow(cosine)

results$results$Cosine = cosine$birth
dif = abs(results$results$Cosine - results$results$Season)
results$results$Difference = apply(cbind(dif,(1-dif)),1,min)
results$results$R = cosine$Pearson
rrr.r = results$results[order(results$results$Cluster),]

# Storing the birth seasonality estimates from two methods
c.rousay = results$results$Cosine
s.rousay = results$results$Season


# ---------------------------------------------------------
# Plot the birth seasonality estimates for the two methods (for above three data)
# ---------------------------------------------------------

dev.off()
par(mar=c(1.5,1.5,1,1), mfrow=c(1,3),oma = c(4, 4, 1, 1))
plot(c.rousay,s.rousay,xlim = c(0.15,0.40),ylim = c(0.15,0.40),pch = 19)
abline(0,1)
legend("topleft","Rousay")
plot(c.bordusani,s.bordusani,xlim = c(0.15,0.65),ylim = c(0.15,0.65),pch = 19)
abline(0,1)
legend("topleft","Bordusani")
plot(c.armenia,s.armenia,xlim = c(0,1),ylim = c(0,1),pch = 19)
abline(0,1)
legend("topleft","Tsaghkahovit")
title(xlab = "estimated birth seasonality from cosine method",ylab = "estimated birth seasonality from SCEM",outer = TRUE, line = 2)


# ---------------------------------------------------------
# Use Armenia data to show the sensitivity to initialization (for Cosine method) 
# ---------------------------------------------------------

# Fix the initial values of amplitude and intercept
amp = seq(1,10,by=0.5)
int = seq(-25,0,by=0.5)

# Find the estimate of birth seasonality, period, delay and the R^2 for all cases
cosfit.birth = array(0,dim = c(length(amp),length(int),length(armenia)))
cosfit.period = array(0,dim = c(length(amp),length(int),length(armenia)))
cosfit.delay = array(0,dim = c(length(amp),length(int),length(armenia)))
cosfit.R = array(0,dim = c(length(amp),length(int),length(armenia)))
for (i in 1:length(amp)){
  for (j in 1:length(int)){
    dd = makeFitsWrong(armenia,amp[i],int[j])
    cosfit.birth[i,j,] = dd$birth
    cosfit.period[i,j,] = dd$X
    cosfit.delay[i,j,] = dd$x0
    cosfit.R[i,j,] = dd$Pearson
  }
}
mm = numeric(length(armenia))
vv = mm
for (i in 1:length(armenia)){
  mm[i] = mean(c(cosfit.birth[,,i]))
  vv[i] = var(c(cosfit.birth[,,i]))
}

# Plots of the model fits to show the sensitivity to initial conditions
dev.off()
par(mfrow = c(2,2),mar = c(4,4,1,1))
plot(mm,sqrt(vv),xlab = "mean of birth season",ylab = "standard deviation",pch = 19)
hist(c(cosfit.period),25,freq = F,xlab = "estimated period",main = "")
hist(c(cosfit.delay),25,freq = F,xlab = "estimated delay",main = "")
hist(c(cosfit.R^2),25,freq = F,xlab = "R-squared value",main = "")


# ---------------------------------------------------------
# Simulate new data (from Rousay data) to compare performance of SCEM and Cosine method
# ---------------------------------------------------------

# Fix the proportion and number of time series to remove data points
removalp = c(0.1,0.2,0.5)
removaln = c(1,3,5,10)

# Find the errors (MSE) and overall variation (MV) for all cases
scam.err = matrix(nrow = 3,ncol = 4)
cosine.err = matrix(nrow = 3,ncol = 4)
scam.var = matrix(nrow = 3,ncol = 4)
cosine.var = matrix(nrow = 3,ncol = 4)
for (i in 1:3){
  for (j in 1:4){
    nn = removaln[j]
    pp = removalp[i]
    temp = matrix(nrow = 100,ncol = 4)
    for (turn in 1:100){
      dd = rousay
      ss = sample(length(dd),nn)
      for (k in 1:nn){
        ii = ceiling(nrow(dd[[ss[k]]])*(1-pp))
        dd[[ss[k]]] = dd[[ss[k]]][1:ii,]
      }
      res = SCEM(dd,-0.33)
      temp[turn,1] = sum((0.28-res$results$Season)^2,na.rm = T)/(length(dd))
      temp[turn,3] = var(res$results$Season)
      cosine = makeFits(dd)
      temp[turn,2] = sum((0.28-cosine$birth)^2,na.rm = T)/(length(dd)-sum(is.na(0.28-cosine$birth)))
      temp[turn,4] = var(cosine$birth)
    }
    scam.err[i,j] = mean(temp[,1])
    cosine.err[i,j] = mean(temp[,2])
    scam.var[i,j] = mean(temp[,3])
    cosine.var[i,j] = mean(temp[,4])
  }
}

