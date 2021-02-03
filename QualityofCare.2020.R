setwd("~/COVID/Quality of Care")
#load libraries
library(deSolve)
library(viridis)

#Set up parameters
parameters.default <- c(gamm1 = 1/2, gamm2 = 1/2, bet1 = 1, bet2 = 1, vac = .02 , pref = .5, rho = .01, mid = .5, sigmoid = 4.5)
parameters = parameters.default
parameters.1p <- parameters.default[-c(2, 4)]

#Set up states, initial values
state <- c(S=9799, I=1, R = 0, V = 0, M = 200, L = 0, H = 0, W = 0, Z = 0, D = 0)
state.1p <- c(S=9999, I=1, R=0, V=0, Z=0, D = 0)

#simplest model - Alt. model 3
alt3 <-function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    #find transitions
    N = S+I+R+V #Population over time
    dS = -bet1 * S * (I)/(N) -  pref*vac * S    #differential equation for Susceptible gen pop Class
    dI = bet1 * S * (I)/(N) - (gamm1) * (I/(1-rho))   #differential equation for Infected gen pop Class
    dR = gamm1 * I         #differential equation for Recovered gen pop Class
    dV = pref*vac * S      #differential equation for Vaccinated gen pop Class
    dZ = bet1 * S *(I)/(N)
    dD = rho * I * gamm1/(1 - rho)
    #return transitions
    list(c(dS, dI, dR, dV, dZ, dD))   
  })
}

#model with Q but not HCW model - Alt. model 1
alt1 <-function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    #find transitions
    N = S+I+R+V #Population over time
    realprop = (S+R+V)/10000 #prop of healthy HCW
    Q = (1/(1+exp(-sigmoid*(realprop-(1-mid))))) #Introducting quality of care
    Qmax = (1/(1+exp(-sigmoid*(mid)))) #Normalizing
    Qmin = (1/(1+exp(-sigmoid*(mid-1)))) #Normalizing
    Q = (Q - Qmin)/(Qmax - Qmin)
    dS = -bet1 * S * (I)/(N) - (Q) * pref * vac * S    #differential equation for Susceptible gen pop Class
    dI = bet1 * S * (I)/(N) - (gamm1) * (I/(1-rho))   #differential equation for Infected gen pop Class
    dR = gamm1 * I * (Q)         #differential equation for Recovered gen pop Class
    dV = pref* Q * vac * S      #differential equation for Vaccinated gen pop Class
    dZ = bet1 * S * (I)/(N)
    dD = (gamm1*I - Q*gamm1*I + Q*rho*gamm1*I)/(1 - rho)
    #return transitions
    list(c(dS, dI, dR, dV, dZ, dD))   
  })
}

#model with HCW separated but not Q - Alt. model 2
alt2 <-function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    #find transitions
    N = S+I+R+V+M+L+H+W #Population over time
    dS = -bet1 * S * (I + L)/(N) -  pref * vac * S    #differential equation for Susceptible gen pop Class
    dI = bet1 * S * (I + L)/(N) - (gamm1) * (I/(1-rho))   #differential equation for Infected gen pop Class
    dR = gamm1 * I          #differential equation for Recovered gen pop Class
    dV = pref *  vac * S      #differential equation for Vaccinated gen pop Class
    dM = -bet2 * M * (I + L)/(N) - vac * M    #differential equation for Susceptible HCW Class
    dL = bet2 * M * (I + L)/(N) - gamm2 * L/(1 - rho)     #differential equation for Infected HCW Class
    dH = gamm2 * L            #differential equation for Recovered HCW Class
    dW =  vac * M             #differential equation for Vaccinated HCW Class
    dZ = bet1 * S * (I + L)/(N) + bet2 * M * (I + L)/(N)
    dD = (gamm1*I - gamm1*I + rho*gamm1*I + gamm2*L - gamm2*L + rho*gamm2*L)/(1 - rho)
    #return transitions
    list(c(dS, dI, dR, dV, dM, dL, dH, dW, dZ, dD))   
  })
}

#Set up function (Equations for full model, HCW sep population with proportional care)
full <-function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    #find transitions
    N = S+I+R+V+M+L+H+W #Population over time
    realprop = (M+H+W)/200 #prop of healthy HCW
    Q = (1/(1+exp(-sigmoid*(realprop-(1-mid))))) #Introducting quality of care
    Qmax = (1/(1+exp(-sigmoid*(mid)))) #Normalizing
    Qmin = (1/(1+exp(-sigmoid*(mid-1)))) #Normalizing
    Q = (Q - Qmin)/(Qmax - Qmin)
    dS = -bet1 * S * (I + L)/(N) - (Q) * pref * vac * S    #differential equation for Susceptible gen pop Class
    dI = bet1 * S * (I + L)/(N) - (gamm1) * (I/(1-rho))   #differential equation for Infected gen pop Class
    dR = gamm1 * I * (Q)         #differential equation for Recovered gen pop Class
    dV = pref * Q * vac * S      #differential equation for Vaccinated gen pop Class
    dM = -bet2 * M * (I + L)/(N) - (Q) * vac * M    #differential equation for Susceptible HCW Class
    dL = bet2 * M * (I + L)/(N) - gamm2 * L/(1 - rho)     #differential equation for Infected HCW Class
    dH = gamm2 * L * Q           #differential equation for Recovered HCW Class
    dW = Q * vac * M             #differential equation for Vaccinated HCW Class
    dZ = bet1 * S * (I + L)/(N) + bet2 * M * (I + L)/(N) 
    dD = (gamm1*I - Q*gamm1*I + rho*Q*gamm1*I + gamm2*L - Q*gamm2*L + rho*Q*gamm2*L)/(1 - rho)
    #return transitions
    list(c(dS, dI, dR, dV, dM, dL, dH, dW, dZ, dD))   
  })
}

#Set up time vector
times <- seq(0, 100, length = 100)

#run solver
out <- ode(y = state, times = times, func = full, parms = parameters)
out<- as.data.frame(out)

##out is a matrix with 9 columns - time, S, I, R, V, M, L, H, W, and 10001 rows - one for each timestep (100 time units, subdivided into hundredths)

#episize[j, k] = 1000 - (out[nrow(out), "S"]+out[nrow(out), "M"]+out[nrow(out), "V"]+out[nrow(out), "W"]) #Calculation of Episize
#mortsize[j, k] = 1000 - (out[nrow(out), "S"]+out[nrow(out), "I"]+out[nrow(out), "R"]+out[nrow(out), "V"]+out[nrow(out), "M"]
 #                        +out[nrow(out), "L"]+out[nrow(out), "H"]+out[nrow(out), "W"]) #Calculation of Mortsize


qOfCare <- function(realprop, sigmoid, mid){
  Q = (1/(1+exp(-sigmoid*(realprop-(1-mid))))) #Introducting quality of care
  Qmax = (1/(1+exp(-sigmoid*(mid)))) #Normalizing
  Qmin = (1/(1+exp(-sigmoid*(mid-1)))) #Normalizing
  Q = (Q - Qmin)/(Qmax - Qmin)
  return(Q)
}
realprop = seq(0, 1, by = .01)

##FigS1
png("FigS1.png", width = 10, units = "in", height = 8, res = 500)
plot(realprop, qOfCare(1-realprop, 4.5, .9), type = "l", lwd = 2, xlab = "Proportion of Infected HCWs", ylab = "Quality of Care")
lines(realprop, qOfCare(1-realprop, 4.5, .1), lwd = 2, lty = 3)
lines(realprop, qOfCare(1-realprop, .5, .9), lwd = 2, col = "grey")
lines(realprop, qOfCare(1-realprop, .5, .1), lwd = 2, col = "grey", lty = 3)
legend("topright", col = c("grey", "grey", "black", "black"), lty = c(1, 3, 3, 1), legend = 
         c("Low Loss Imp (k = 0.5) / High Red (m = 0.9)", "Low Loss Imp (k = 0.5) / Low Red (m = 0.1)", 
           "High Loss Imp (k = 4.5) / Low Red (m = 0.1)", "High Loss Imp (k = 4.5) / High Red (m = 0.9)"), 
       bty = "n", lwd = 2)
dev.off()

#colors
cols = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
##FigS2
png("FigS2.png", width = 10, height = 10,units = "in", res = 500)
layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))
plot(realprop, qOfCare(1-realprop, 4.5, 0.9), type = "l", lwd = 2, xlab = "Proportion of Infected HCWs", 
     ylab = "Quality of Care", col = cols[1])
lines(realprop, qOfCare(1-realprop, 3, 0.9), lwd = 2, lty = 2, col = cols[2])
lines(realprop, qOfCare(1-realprop, 2, 0.9), lwd = 2, col = cols[3], lty = 3)
lines(realprop, qOfCare(1-realprop, 0.5, 0.9), lwd = 2, col = cols[4], lty = 4)
legend("topright", col = cols, lty = 1:4, bty = "n", lwd = 2,
       legend = c("High Loss Impact (k = 4.5)", "Medium-high Loss Impact (k = 3)", 
                  "Medium-low Loss Impact (k = 2)", "Low Loss Impact (k = 0.5)"))
mtext("A", side = 2, line = 1, at = 1.15, cex = 2, family= "serif", las = 1)

plot(realprop, qOfCare(1-realprop, 4.5, .1), type = "l", lwd = 2, xlab = "Proportion of Infected HCWs", 
     ylab = "Quality of Care", col = cols[1])
lines(realprop, qOfCare(1-realprop, 3, .1), lwd = 2, lty = 2, col = cols[2])
lines(realprop, qOfCare(1-realprop, 2, .1), lwd = 2, col = cols[3], lty = 3)
lines(realprop, qOfCare(1-realprop, .5, .1), lwd = 2, col = cols[4], lty = 4)
legend("topright", col = cols, lty = 1:4, bty = "n", lwd = 2,
       legend = c("High Loss Impact (k = 4.5)", "Medium-high Loss Impact (k=3)", 
                  "Medium-low Loss Impact (k = 2)", "Low Loss Impact (k = 0.5)"))
mtext("B", side = 2, line = 1, at = 1.15, cex = 2, family= "serif", las = 1)

plot(realprop, qOfCare(1-realprop, 4.5, .9), type = "l", lwd = 2, xlab = "Proportion of Infected HCWs", 
     ylab = "Quality of Care", col = cols[1])
lines(realprop, qOfCare(1-realprop, 4.5, .7), lwd = 2, lty = 2, col = cols[2])
lines(realprop, qOfCare(1-realprop, 4.5, .3), lwd = 2, col = cols[3], lty = 3)
lines(realprop, qOfCare(1-realprop, 4.5, .1), lwd = 2, col = cols[4], lty = 4)
legend("topright", col = cols, lty = 1:4, bty = "n", lwd = 2,
       legend = c("High Redundancy (m = 0.9)", "Medium-high Redundancy (m = 0.7)", 
                  "Medium-low Redundancy (m = 0.3)", "Low Redundancy (m = 0.1)"))
mtext("C", side = 2, line = 1, at = 1.15, cex = 2, family= "serif", las = 1)

plot(realprop, qOfCare(1-realprop, .5, .9), type = "l", lwd = 2, xlab = "Proportion of Infected HCWs", 
     ylab = "Quality of Care", col = cols[1])
lines(realprop, qOfCare(1-realprop, .5, .7), lwd = 2, lty = 2, col = cols[2])
lines(realprop, qOfCare(1-realprop, .5, .3), lwd = 2, col = cols[3], lty = 3)
lines(realprop, qOfCare(1-realprop, .5, .1), lwd = 2, col = cols[4], lty = 4)
legend("topright", col = cols, lty = 1:4, bty = "n", lwd = 2,
       legend = c("High Redundancy (m = 0.9)", "Medium-high Redundancy (m = 0.7)", 
                  "Medium-low Redundancy (m = 0.3)", "Low Redundancy (m = 0.1)"))
mtext("D", side = 2, line = 1, at = 1.15, cex = 2, family= "serif", las = 1)
dev.off()

##Fig 2
times <- seq(0, 100, length = 100)
png("Fig2.png", width = 10, height = 5,units = "in", res = 500)
#run solver
layout(matrix(1:2, nrow = 1))
parameters = parameters.default
out <- ode(y = state, times = times, func = alt2, parms = parameters)
out1<- as.data.frame(out)
parameters = parameters.default
parameters[names(parameters) == "mid"] = 0
out <- ode(y = state, times = times, func = full, parms = parameters)
out2<- as.data.frame(out)
parameters[names(parameters) == "mid"] = 1
out <- ode(y = state, times = times, func = full, parms = parameters)
out3<- as.data.frame(out)
plot(NA, xlim = c(0, 100), ylim = c(0, 7000), xlab = "Weeks", ylab = "Infected")
polygon(c(rev(out2$time), out3$time), c(rev(out2$Z), out3$Z), col = "grey", border = NA)
lines(out1$time, out1$Z, lwd = 2)
legend("bottomright", lwd = c(2, 4), col = c("black", "grey"), bty = "n", legend = 
         c("No Proportional Care", "Proportional Care"))
mtext("A", side = 2, line = 1, at = 8000, family = "serif", las = 1, cex = 2)

parameters = parameters.default
parameters[names(parameters)]
parameters[names(parameters) == "mid"] = 1
parameters[names(parameters) == "bet2"] = 1.5
outfull <- ode(y = state, times = times, func = full, parms = parameters)
out.full.1 <- as.data.frame(outfull)
parameters.1p <- parameters[-c(2, 4)]
outalt1 <- ode(y = state.1p, times = times, func = alt1, parms = parameters.1p)
out.alt1.1 <- as.data.frame(outalt1)
parameters[names(parameters) == "mid"] = 0
outfull <- ode(y = state, times = times, func = full, parms = parameters)
out.full.2 <- as.data.frame(outfull)
parameters.1p <- parameters[-c(2, 4)]
outalt1 <- ode(y = state.1p, times = times, func = alt1, parms = parameters.1p)
out.alt1.2 <- as.data.frame(outalt1)
outalt2 <- ode(y = state, times = times, func = alt2, parms = parameters)
out.alt2 <- as.data.frame(outalt2)
outalt3 <- ode(y = state.1p, times = times, func = alt3, parms = parameters.1p)
out.alt3 <- as.data.frame(outalt3)
plot(NA, xlim = c(0, 100), ylim = c(0, 7000), xlab = "Weeks", ylab = "Infected")
polygon(c(rev(out.full.1$time), out.full.2$time), c(rev(out.full.1$Z), out.full.2$Z), col = "grey80", border = NA)
polygon(c(rev(out.alt1.1$time), out.alt1.2$time), c(rev(out.alt1.1$Z), out.alt1.2$Z), col = "grey50", border = NA)
lines(out.alt2$time, out.alt2$Z, col = "grey20", lwd = 2)
lines(out.alt3$time, out.alt3$Z, col = "black", lwd = 2, lty = 2)
legend("bottomright", col = c("grey80", "grey50", "grey20", "black"), lwd = c(4, 4, 2, 2), lty = c(1, 1, 1, 2), 
       bty = "n", legend = c("Full Model","Alternate Model I", "Alternate Model II", "Alternate Model III"))
mtext("B", side = 2, line = 1, at = 8000, family = "serif", las = 1, cex = 2)
dev.off()

###### Everything below this involves looping through different values of midpoint and sigmoid parameter
midpointparameter = seq(0, 1, length = 100) #values of midpoint sigmoid parameter
sigmoidparameter = seq(.1, 5, length = 100) #values of sigmoid slope parameter

episize = matrix(rep(0, 100*100), ncol = 100, nrow = 100) #Creates a matrix for the epidemic size thats 100 by 100
mortsize = matrix(rep(0, 100*100), ncol = 100, nrow = 100) #Creates a matrix for the mortality size thats 100 by 100

for(j in 1:100){ #loop iterating all possible midpoint rates
  for(k in 1:100){ #loop iterating all possible sigmoidal/slope parameters
    #Set up parameters
    parameters = parameters.default
    parameters[names(parameters) == "mid"] = midpointparameter[j]
    parameters[names(parameters) == "sigmoid"] = sigmoidparameter[k]
    out <- ode(y = state, times = times, func = full, parms = parameters)
    out<- as.data.frame(out)
    
    ##out is a matrix with 9 columns - time, S, I, R, V, M, L, H, W, and 10001 rows - one for each timestep (100 time units, subdivided into hundredths)
    
    episize[j, k] = out$Z[nrow(out)]
    mortsize[j, k] = out$D[nrow(out)]  
  } #end of loop of k (slope/sigmoidal parameter)
} #end of loop of j (midpoint parameter)

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='', ax = TRUE) {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(0.5,4), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  if(ax){
    axis(2, ticks, las=1, cex.axis=1)
  }
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,4,y+1/scale, col=lut[i], border=NA)
  }
}

#Plot dimensions and margins
png("Fig3.png", width = 16, height = 5,units= "in", res = 500)
layout(matrix(c(1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 5, 6), nrow = 1))
par(mar=c(6,5,5,3))

image(x=midpointparameter, y = sigmoidparameter, z = log(1+episize), breaks = log(1+seq(5000, 6000, length = 101)), 
      col = rev(inferno(100)), xlab = "Redundancy Parameter (m)", ylab = "Loss Impact Parameter (k)",
      cex.lab = 1.5, cex.main = 3) #creates image of episize
contour(x=midpointparameter, y = sigmoidparameter, 
        levels = seq(0, 10000, by = 50),  
        z = episize, add = TRUE, col='black') #adds contour lines to the image, showing all combinations of recovery rates that lead to a distinct episize by the hundreds from 100 to 900
abline(v = 0.5, col = "white", lty = 2)
abline(h = 2.5, col = "white", lty = 2)
mtext("A", side = 2, line = 1, at = 5.5, family = "serif", las = 1, cex = 2)
color.bar(lut = rev(inferno(100)), min = log(1+5000), max = log(1+6000), title = "Infected", ax = FALSE)
axis(2, at = log(1+ seq(5000, 6000, length = 11)), labels = seq(5000, 6000, length = 11), las = 1)

image(x=midpointparameter, y = sigmoidparameter, z = mortsize, col = rev(inferno(100)), 
      breaks = seq(70, 1600, length = 101), xlab = "Redundancy Parameter (m)", ylab = "Loss Impact Parameter (k)", 
      cex.lab = 1.5, cex.main = 3)
contour(x=midpointparameter, y = sigmoidparameter, z =  mortsize, add = TRUE, 
        levels = seq(50, 1600, by = 50), col='black')
abline(v = 0.5, col = "white", lty = 2)
abline(h = 2.5, col = "white", lty = 2)
mtext("B", side = 2, line = 1, at = 5.5, family = "serif", las = 1, cex = 2)
color.bar(lut = rev(inferno(100)), min = 70, max = 1600, title = "Dead")

cfr = mortsize/episize
cfr[which(episize == 0)[1], which(episize == 0)[2]] = 0
image(x=midpointparameter, y = sigmoidparameter, z = cfr, col = rev(inferno(100)),
      xlab = "Redundancy Parameter (m)", ylab = "Loss Impact Parameter (k)", breaks = seq(0, .26, length = 101),
      cex.lab = 1.5, cex.main = 3)
contour(x=midpointparameter, y = sigmoidparameter, z = 100*cfr, add = TRUE, levels = seq(0, 26, by = 1), col='black')
abline(v = 0.5, col = "white", lty = 2)
abline(h = 2.5, col = "white", lty = 2)
mtext("C", side = 2, line = 1, at = 5.5, family = "serif", las = 1, cex = 2)
color.bar(lut = rev(inferno(100)), min = 0, max = 26, title = "Percent\nMortality", nticks = 14)
dev.off()
###

#Set up function (Equations for full model, HCW sep population with proportional care)
SIR2pq1 <-function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    #find transitions
    N = S+I+R+V+M+L+H+W #Population over time
    realprop = (M+H+W)/200 #prop of healthy HCW
    Q = 1
    dS = -bet1 * S * (I + L)/(N) - (Q) * pref * vac * S    #differential equation for Susceptible gen pop Class
    dI = bet1 * S * (I + L)/(N) - (gamm1) * (I/(1-rho))   #differential equation for Infected gen pop Class
    dR = gamm1 * I * (Q)         #differential equation for Recovered gen pop Class
    dV = pref * Q * vac * S      #differential equation for Vaccinated gen pop Class
    dM = -bet2 * M * (I + L)/(N) - (Q) * vac * M    #differential equation for Susceptible HCW Class
    dL = bet2 * M * (I + L)/(N) - gamm2 * L/(1 - rho)     #differential equation for Infected HCW Class
    dH = gamm2 * L * Q           #differential equation for Recovered HCW Class
    dW = Q * vac * M             #differential equation for Vaccinated HCW Class
    dZ = bet1 * S * (I + L)/(N) + bet2 * M * (I + L)/(N) 
    dD = rho * (I * gamm1 + L * gamm2)/(1 - rho)
    #return transitions
    list(c(dS, dI, dR, dV, dM, dL, dH, dW, dZ, dD))   
  })
}
png("Fig4.png", width = 16, height = 8, units = "in", res = 500)
cols = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
layout(matrix(1:2, nrow = 1))
parameters <- parameters.default
parameters[names(parameters) == "mid"] = .1
parameters[names(parameters) == "sigmoid"] = 4.5
out <- ode(y = state, times = times, func = full, parms = parameters)
out<- as.data.frame(out)
plot(out$time, (out$M + out$H + out$W)/200, type = "l", lwd = 2, xlab = "Weeks", 
     ylab = "Proportion of Population", ylim = c(0, 1), col=cols[1])
lines(out$time, (out$V + out$W)/10000, lty = 2, col = cols[2], lwd = 2)
lines(out$time, (out$D)/10000, lty = 2, col = cols[3], lwd = 2)
lines(out$time, (out$R + out$H)/10000, lty = 2, col = cols[4], lwd = 2)
out2 <- ode(y = state, times = times, func = SIR2pq1, parms = parameters)
out2<- as.data.frame(out2)
lines(out2$time, (out2$V + out2$W)/10000, col = cols[2], lty = 3, lwd = 2)
lines(out2$time, (out2$D)/10000, lty = 3, col = cols[3], lwd = 2) #"#41ab5d"
lines(out2$time, (out2$R + out2$H)/10000, lty = 3, col = cols[4], lwd = 2)
legend(50, .85, col = c(cols, cols[2], cols[3],cols[4]), lwd = 2, lty = c(1, 2, 2, 2, 3, 3, 3), bty = "n",
       legend = c("Healthy HCW", "Cumulative Vaccinated", "Cumulative Deaths", "Cumulative Recovered", 
                  "Cumulative Vaccinated (Q = 1)", "Cumulative Deaths (Q = 1)", "Cumulative Recovered (Q = 1)"))
mtext("A", side = 2, line = 1, at = 1.15, family = "serif", las = 1, cex = 2)

plot(NA, xlim = c(0, 100), ylim = c(0, .2), xlab = "Weeks", ylab = "Proportion of the Population")
out <- ode(y = state, times = times, func = SIR2pq1, parms = parameters)
out<- as.data.frame(out)
lines(out$time, out$D/10000, lwd = 2, col = "grey75", lty = 3)
parameters[names(parameters) == "mid"] = .9
parameters[names(parameters) == "sigmoid"] = 4.5
out <- ode(y = state, times = times, func = full, parms = parameters)
out<- as.data.frame(out)
lines(out$time, out$D/10000, lwd = 2, col = cols[1], lty = 3)
parameters[names(parameters) == "mid"] = .9
parameters[names(parameters) == "sigmoid"] = .5
out <- ode(y = state, times = times, func = full, parms = parameters)
out<- as.data.frame(out)
lines(out$time, out$D/10000, lwd = 2, col = cols[2], lty = 3)
parameters[names(parameters) == "mid"] = .1
parameters[names(parameters) == "sigmoid"] = 4.5
out <- ode(y = state, times = times, func = full, parms = parameters)
out<- as.data.frame(out)
lines(out$time, out$D/10000, lwd = 2, col = cols[3], lty = 3)
parameters[names(parameters) == "mid"] = .1
parameters[names(parameters) == "sigmoid"] = .5
out <- ode(y = state, times = times, func = full, parms = parameters)
out<- as.data.frame(out)
lines(out$time, out$D/10000, lwd = 2, col = cols[4], lty = 3)
legend("topright", bty = "n", lty = 3, lwd = 2, col = c("grey75", cols), 
       legend = c("Cumulative Deaths, Q = 1", "Cumulative Deaths,\nHigh Loss Impact, High Redundancy",
                  "Cumulative Deaths,\nLow Loss Impact, High Redundancy", "Cumulative Deaths,\nHigh Loss Impact, Low Redundancy", 
                  "Cumulative Deaths,\nLow Loss Impact, Low Redundancy"), y.intersp = 2)
mtext("B", side = 2, line = 1, at = .23, family = "serif", las = 1, cex = 2)
dev.off()


parameters.default <- c(gamm1 = 1/2, gamm2 = 1/2, bet1 = 1, bet2 = 1, vac = .02 , pref = .5, rho = .01, mid = .9, sigmoid = .5)
parameters = parameters.default
state <- c(S=9799, I=1, R = 0, V = 0, M = 200, L = 0, H = 0, W = 0, Z = 0, D = 0)

death <- function(parameters, state){
  out <- ode(y = state, times = times, func = full, parms = parameters)
  out<- as.data.frame(out)
  return(out$D[length(out$D)])  
}
cases <- function(parameters, state){
  out <- ode(y = state, times = times, func = full, parms = parameters)
  out<- as.data.frame(out)
  return(out$Z[length(out$Z)])  
}
deaths0.9.ll.hr <- rep(NA, length(parameters))
deaths1.1.ll.hr <- rep(NA, length(parameters))
deaths0.9.ll.lr <- rep(NA, length(parameters))
deaths1.1.ll.lr <- rep(NA, length(parameters))
deaths0.9.hl.hr <- rep(NA, length(parameters))
deaths1.1.hl.hr <- rep(NA, length(parameters))
deaths0.9.hl.lr <- rep(NA, length(parameters))
deaths1.1.hl.lr <- rep(NA, length(parameters))

cases0.9.ll.hr <- rep(NA, length(parameters))
cases1.1.ll.hr <- rep(NA, length(parameters))
cases0.9.ll.lr <- rep(NA, length(parameters))
cases1.1.ll.lr <- rep(NA, length(parameters))
cases0.9.hl.hr <- rep(NA, length(parameters))
cases1.1.hl.hr <- rep(NA, length(parameters))
cases0.9.hl.lr <- rep(NA, length(parameters))
cases1.1.hl.lr <- rep(NA, length(parameters))

for(i in 1:7){
  mod <- rep(1, length(parameters))
  mod[i] <- 0.9
  parameters = parameters.default
  deaths0.9.ll.hr[i] <- death(mod*parameters, state)
  cases0.9.ll.hr[i] <- cases(mod*parameters, state)
  mod[i] <- 1.1
  deaths1.1.ll.hr[i] <- death(mod*parameters, state)
  cases1.1.ll.hr[i] <- cases(mod*parameters, state)
  mod[i] <- 0.9
  parameters = parameters.default
  parameters[8] = 0.1
  deaths0.9.ll.lr[i] <- death(mod*parameters, state)
  cases0.9.ll.lr[i] <- cases(mod*parameters, state)
  mod[i] <- 1.1
  deaths1.1.ll.lr[i] <- death(mod*parameters, state)
  cases1.1.ll.lr[i] <- cases(mod*parameters, state)
  mod[i] <- 0.9
  parameters = parameters.default
  parameters[9] = 4.5
  deaths0.9.hl.hr[i] <- death(mod*parameters, state)
  cases0.9.hl.hr[i] <- cases(mod*parameters, state)
  mod[i] <- 1.1
  deaths1.1.hl.hr[i] <- death(mod*parameters, state)
  cases1.1.hl.hr[i] <- cases(mod*parameters, state)
  mod[i] <- 0.9
  parameters = parameters.default
  parameters[8] = 0.1
  parameters[9] = 4.5
  deaths0.9.hl.lr[i] <- death(mod*parameters, state)
  cases0.9.hl.lr[i] <- cases(mod*parameters, state)
  mod[i] <- 1.1
  deaths1.1.hl.lr[i] <- death(mod*parameters, state)
  cases1.1.hl.lr[i] <- cases(mod*parameters, state)
}

png("death.sens.png", width = 10, height = 10, units = "in", res = 500)
par(
  mfrow=c(2,2), 
  oma = c(0, 0, 2, 0),
  mar=c(4,4,2,2))
  plot(NA, xlim = c(0, 1600), ylim = c(0, 8), main = "Low loss-impact, high redundancy", xlab = "deaths", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(deaths0.9.ll.hr[i] > deaths1.1.ll.hr[i]){
      lines(c(deaths0.9.ll.hr[i], deaths1.1.ll.hr[i]), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(deaths0.9.ll.hr[i] <= deaths1.1.ll.hr[i]){
      lines(c(deaths0.9.ll.hr[i], deaths1.1.ll.hr[i]), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .9
  parameters[9] = .5
  abline(v = death(parameters, state), col = "grey")
  legend("topright", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)

    plot(NA, xlim = c(0, 1600), ylim = c(0, 8), main = "Low loss-impact, low redundancy", xlab = "deaths", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(deaths0.9.ll.lr[i] > deaths1.1.ll.lr[i]){
      lines(c(deaths0.9.ll.lr[i], deaths1.1.ll.lr[i]), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(deaths0.9.ll.lr[i] <= deaths1.1.ll.lr[i]){
      lines(c(deaths0.9.ll.lr[i], deaths1.1.ll.lr[i]), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .1
  parameters[9] = .5
  abline(v = death(parameters, state), col = "grey")
  legend("topright", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)

    plot(NA, xlim = c(0, 1600), ylim = c(0, 8), main = "High loss-impact, high redundancy", xlab = "deaths", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(deaths0.9.hl.hr[i] > deaths1.1.hl.hr[i]){
      lines(c(deaths0.9.hl.hr[i], deaths1.1.hl.hr[i]), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(deaths0.9.hl.hr[i] <= deaths1.1.hl.hr[i]){
      lines(c(deaths0.9.hl.hr[i], deaths1.1.hl.hr[i]), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .9
  parameters[9] = 4.5
  abline(v = death(parameters, state), col = "grey")
  legend("topright", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  plot(NA, xlim = c(0, 1600), ylim = c(0, 8), main = "High loss-impact, low redundancy", xlab = "deaths", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(deaths0.9.hl.lr[i] > deaths1.1.hl.lr[i]){
      lines(c(deaths0.9.hl.lr[i], deaths1.1.hl.lr[i]), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(deaths0.9.hl.lr[i] <= deaths1.1.hl.lr[i]){
      lines(c(deaths0.9.hl.lr[i], deaths1.1.hl.lr[i]), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .1
  parameters[9] = 4.5
  abline(v = death(parameters, state), col = "grey")
  legend("topright", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  
  mtext("Sensitivity of death count in all four regimes", outer = TRUE, cex = 1.5)
  dev.off()
  
  png("cases.sens.png", width = 10, height = 10, units = "in", res = 500)
  par(
    mfrow=c(2,2), 
    oma = c(0, 0, 2, 0),
    mar=c(4,4,2,2))
  plot(NA, xlim = c(0, 7000), ylim = c(0, 8), main = "Low loss-impact, high redundancy", xlab = "cases", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(cases0.9.ll.hr[i] > cases1.1.ll.hr[i]){
      lines(c(cases0.9.ll.hr[i], cases1.1.ll.hr[i]), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(cases0.9.ll.hr[i] <= cases1.1.ll.hr[i]){
      lines(c(cases0.9.ll.hr[i], cases1.1.ll.hr[i]), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .9
  parameters[9] = .5
  abline(v = cases(parameters, state), col = "grey")
  legend("topleft", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  
  plot(NA, xlim = c(0, 7000), ylim = c(0, 8), main = "Low loss-impact, low redundancy", xlab = "cases", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(cases0.9.ll.lr[i] > cases1.1.ll.lr[i]){
      lines(c(cases0.9.ll.lr[i], cases1.1.ll.lr[i]), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(cases0.9.ll.lr[i] <= cases1.1.ll.lr[i]){
      lines(c(cases0.9.ll.lr[i], cases1.1.ll.lr[i]), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .1
  parameters[9] = .5
  abline(v = cases(parameters, state), col = "grey")
  legend("topleft", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  
  plot(NA, xlim = c(0, 7000), ylim = c(0, 8), main = "High loss-impact, high redundancy", xlab = "cases", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(cases0.9.hl.hr[i] > cases1.1.hl.hr[i]){
      lines(c(cases0.9.hl.hr[i], cases1.1.hl.hr[i]), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(cases0.9.hl.hr[i] <= cases1.1.hl.hr[i]){
      lines(c(cases0.9.hl.hr[i], cases1.1.hl.hr[i]), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .9
  parameters[9] = 4.5
  abline(v = cases(parameters, state), col = "grey")
  legend("topleft", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  plot(NA, xlim = c(0, 7000), ylim = c(0, 8), main = "High loss-impact, low redundancy", xlab = "cases", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(cases0.9.hl.lr[i] > cases1.1.hl.lr[i]){
      lines(c(cases0.9.hl.lr[i], cases1.1.hl.lr[i]), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(cases0.9.hl.lr[i] <= cases1.1.hl.lr[i]){
      lines(c(cases0.9.hl.lr[i], cases1.1.hl.lr[i]), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .1
  parameters[9] = 4.5
  abline(v = cases(parameters, state), col = "grey")
  legend("topleft", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  
  mtext("Sensitivity of case count in all four regimes", outer = TRUE, cex = 1.5)
  dev.off()
  
  png("cfr.sens.png", width = 10, height = 10, units = "in", res = 500)
  par(
    mfrow=c(2,2), 
    oma = c(0, 0, 2, 0),
    mar=c(4,4,2,2))
  plot(NA, xlim = c(0, 30), ylim = c(0, 8), main = "Low loss-impact, high redundancy", xlab = "cfr (%)", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(deaths0.9.ll.hr[i]/cases0.9.ll.hr[i] > deaths1.1.ll.hr[i]/cases1.1.ll.hr[i]){
      lines(c(deaths0.9.ll.hr[i]/cases0.9.ll.hr[i]*100, deaths1.1.ll.hr[i]/cases1.1.ll.hr[i]*100), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(deaths0.9.ll.hr[i]/cases0.9.ll.hr[i] <= deaths1.1.ll.hr[i]/cases1.1.ll.hr[i]){
      lines(c(deaths0.9.ll.hr[i]/cases0.9.ll.hr[i]*100, deaths1.1.ll.hr[i]/cases1.1.ll.hr[i]*100), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .9
  parameters[9] = .5
  abline(v = death(parameters, state)/cases(parameters, state)*100, col = "grey")
  legend("topright", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  
  plot(NA, xlim = c(0, 30), ylim = c(0, 8), main = "Low loss-impact, low redundancy", xlab = "cfr (%)", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(deaths0.9.ll.lr[i]/cases0.9.ll.lr[i] > deaths1.1.ll.lr[i]/cases1.1.ll.lr[i]){
      lines(c(deaths0.9.ll.lr[i]/cases0.9.ll.lr[i]*100, deaths1.1.ll.lr[i]/cases1.1.ll.lr[i]*100), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(deaths0.9.ll.lr[i]/cases0.9.ll.lr[i] <= deaths1.1.ll.lr[i]/cases1.1.ll.lr[i]){
      lines(c(deaths0.9.ll.lr[i]/cases0.9.ll.lr[i]*100, deaths1.1.ll.lr[i]/cases1.1.ll.lr[i]*100), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .1
  parameters[9] = .5
  abline(v = death(parameters, state)/cases(parameters, state)*100, col = "grey")
  legend("topright", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  
  plot(NA, xlim = c(0, 30), ylim = c(0, 8), main = "High loss-impact, high redundancy", xlab = "cfr (%)", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(deaths0.9.hl.hr[i]/cases0.9.hl.hr[i] > deaths1.1.hl.hr[i]/cases1.1.hl.hr[i]){
      lines(c(deaths0.9.hl.hr[i]/cases0.9.hl.hr[i]*100, deaths1.1.hl.hr[i]/cases1.1.hl.hr[i]*100), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(deaths0.9.hl.hr[i]/cases0.9.hl.hr[i] <= deaths1.1.hl.hr[i]/cases1.1.hl.hr[i]){
      lines(c(deaths0.9.hl.hr[i]/cases0.9.hl.hr[i]*100, deaths1.1.hl.hr[i]/cases1.1.hl.hr[i]*100), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .9
  parameters[9] = 4.5
  abline(v = death(parameters, state)/cases(parameters, state)*100, col = "grey")
  legend("topright", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
 
   plot(NA, xlim = c(0, 30), ylim = c(0, 8), main = "High loss-impact, low redundancy", xlab = "cfr (%)", yaxt = "n", ylab = "")
  axis(2, at = 1:7, labels = names(parameters)[1:7], las = 2)
  for(i in 1:7){
    if(deaths0.9.hl.lr[i]/cases0.9.hl.lr[i] > deaths1.1.hl.lr[i]/cases1.1.hl.lr[i]){
      lines(c(deaths0.9.hl.lr[i]/cases0.9.hl.lr[i]*100, deaths1.1.hl.lr[i]/cases1.1.hl.lr[i]*100), c(i, i), lwd = 2, col = "#d95f02")
    }
    if(deaths0.9.hl.lr[i]/cases0.9.hl.lr[i] <= deaths1.1.hl.lr[i]/cases1.1.hl.lr[i]){
      lines(c(deaths0.9.hl.lr[i]/cases0.9.hl.lr[i]*100, deaths1.1.hl.lr[i]/cases1.1.hl.lr[i]*100), c(i, i), lwd = 2, col = "#7570b3")
    }
  }
  parameters= parameters.default
  parameters[8] = .1
  parameters[9] = 4.5
  abline(v = death(parameters, state)/cases(parameters, state)*100, col = "grey")
  legend("topright", bty = "n", col = c("#d95f02", "#7570b3"), legend = c("-10% > +10%", "-10% <= +10%"), lwd = 2)
  
  mtext("Sensitivity of case fatality ratio in all four regimes", outer = TRUE, cex = 1.5)
  dev.off()
  