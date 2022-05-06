#Code to accompany Hudson and Creanza 2022, "Modeling how population size drives the evolution of birdsong, a functional cultural trait" https://doi.org/10.1111/evo.14489


# Set up ------------------------------------------------------------------


library(pacman)
p_load(evd,sfsmisc,ggplot2,cowplot,dplyr, doParallel,igraph, pheatmap, RColorBrewer, grid, lattice)
Euler = 0.5772156649

#If you're running all of the script at once, set this to true so figures don't plot on top of each other
new_windows <- F


# Three simulation functions  ---------------------------------------------

#Three different simulation functions: 
#Z.N: Get Z in a small network with no degree/k.
#Use this function to substract from each col in panel 3a. 
simulateZ.N<- function(alph,bet,gens,smallgraf,Npop,...){
  smallgraf <- make_full_graph(n=Npop,directed = F,loops = F)
  
  V(smallgraf)$newZ <- 1 #initialize at 0 for each set of alph and bet
  V(smallgraf)$oldZ <- 1
  Npop=max(V(smallgraf))
  for(i in 1:gens){
    for(j in 1:Npop) {
      
      V(smallgraf)[j]$oldZ<-V(smallgraf)[j]$newZ
      zmax<-max(V(smallgraf)$oldZ)
      newZ<-rgumbel(1,loc=zmax - alph,scale=bet)
      V(smallgraf)[j]$newZ<-newZ
    }#End iterating across all indivs
    ##print(paste("gen",i,"done"))
  }
  return(data.frame(alph=alph,bet=bet,Ngens=gens,meanZ=mean(V(smallgraf)$newZ)))
  
}

#Z.D: Simulate a subset of K (or D for degree) tutors (static throughout all timesteps)
simulateZ.D<- function(alph,bet,gens,Npop,degnet,...){
  V(degnet)$newZ <- 1 #initialize at 1 for each set of alph and bet
  V(degnet)$oldZ <- 1
  for(i in 1:gens){
    for(currIndiv in 1:Npop)
    {
      V(degnet)[currIndiv]$oldZ<-V(degnet)[currIndiv]$newZ
      connectednodes <- neighbors(degnet, currIndiv)
      subs_zmax <- max(V(degnet)$oldZ[connectednodes])
      newZ<-rgumbel(1,loc=subs_zmax - alph,scale=bet)
      V(degnet)[currIndiv]$newZ<-newZ
    }#End iterating across all indivs
  }
  #return(degnet)
  return(data.frame(alph=alph,bet=bet,Ngens=gens,meanZ=mean(V(degnet)$newZ)))
  ##print(paste("gen",i,"mean z",meanZ))
}


#Z.K: Sample K connections out of N (every generation new)
simulateZ.K<- function(alph,bet,gens,K,fullGraf){
  V(fullGraf)$newZ <- 1 #initialize at 0 for each set of alph and bet
  V(fullGraf)$oldZ <- 1
  for(i in 1:gens){
    Npop = max(V(fullGraf))
    for(currIndiv in 1:Npop)
    {
      V(fullGraf)[currIndiv]$oldZ<-V(fullGraf)[currIndiv]$newZ
      
      k_indices<-sample(1:vcount(fullGraf),K) #each generation, sample new tutors (k)
      subs_zmax<-max(V(fullGraf)$oldZ[k_indices]) #find max z in that subset
      newZ<-rgumbel(1,loc=subs_zmax - alph,scale=bet)
      V(fullGraf)[currIndiv]$newZ<-newZ
    }#End iterating across all indivs
    #print(paste("Pop avg z=",mean(V(fullGraf)$newZ)))
  }
  funresults <- return(data.frame(alph=alph,bet=bet,Ngens=gens,meanZ=mean(V(fullGraf)$newZ)))
}


# Alpha bottleneck simulation function ---------------------------------------------

AlphaBottleneckSim <-
  function(StartZ,
           alpha1,
           beta1,
           simreps,
           phase1gens = 100,
           phase2gens = 100,
           phase3gens = 0,
           phase1Pop = 1000,
           phase2Pop = 100,
           phase3Pop = 0,
           alpharange = 2,
           verbose = F) {
    BottlePops <-
      c(
        rep(phase1Pop, phase1gens),
        rep(phase2Pop, phase2gens),
        rep(phase3Pop, phase3gens)
      ) 
    simgens <- length(BottlePops)
    alphatable1 <- seq(-1,alpha1*2,length.out = alpha1*20+1)
    squares <-lapply(seq(1,alpha1*5), function (x) x^2)
    squaresToAdd <- as.numeric(c(rep(0,alpha1*15+1),squares))
    alphatable <- as.numeric(alphatable1)+squaresToAdd
    alphatable
    
    Zresults <- rep(0, simgens)
    simmat <-
      matrix(0, nrow = simgens, ncol = simreps)
    Newdist <- 0
    
    Gdist <-
      rgumbel(phase1Pop, loc = (StartZ - alpha1), scale = beta1) #Initial probability of matching StartZ
    alphamat <- matrix(0, nrow = simgens, ncol = simreps)
    maxmat <- matrix(0, nrow = simgens, ncol = simreps)
    indexmat <- matrix(0, nrow = simgens, ncol = simreps)
    
    for (i in 1:simreps) {
      #starting new replicate run at gen 0, re-initialize all values
      Newdist <- Gdist
      currAlpha <- alpha1
      newAlpha <- alpha1
      NewMaxZ <- StartZ
      for (j in 1:simgens) {
        currPop <- BottlePops[j]
        currAlpha <- newAlpha
        Newdist <-
          rgumbel(currPop,
                  loc = (NewMaxZ - currAlpha),
                  scale = beta1)
        NewMaxZ <- max(Newdist)
        
        NewmeanZ <- mean(Newdist)
        Zresults[j] <- NewMaxZ
        simmat[j, i] <- NewmeanZ
        if (NewMaxZ < 0) { 
          NewMaxZ <- 0
        }
        if (NewMaxZ - StartZ <= 0) {
         
          alphaindex <-
            round(length(alphatable) / 2 + (NewMaxZ - StartZ))
          
          #subtract
          if (verbose == T) {
            print(
              paste(
                "NewMaxZ minus start Z: ",
                NewMaxZ - StartZ,
                ", alpha index:",
                alphaindex
              )
            )
          }
          #end of "decrease alpha" condition
        }else{
          #If new Z is LARGER than starting Z (trait is harder) then alpha should increase (larger alpha index than 50)
          #term is postive
          alphaindex <-
            round(length(alphatable) / 2 + (NewMaxZ - StartZ))
          if (verbose == T) {
            print(
              paste(
                "NewMaxZ minus start Z: ",
                NewMaxZ - StartZ,
                ", alpha index:",
                alphaindex
              )
            )
          }
        }#End of "increase alpha" condition
        
        if (alphaindex > length(alphatable)) {
          alphaindex <- length(alphatable)
        }
        if(alphaindex == 0){
          alphaindex <- 1
        }
        
        
        newAlpha <- alphatable[alphaindex]
        indexmat[j, i] <- alphaindex
        alphamat[j, i] <-
          newAlpha
       #print(paste0("New alpha=", newAlpha))
        
        maxmat[j, i] <- NewMaxZ
       
        if (verbose == T) {
          print(paste("gen ", j, ", newAlpha =", newAlpha))
        }
      }   #End of each generation loop
      
      if (verbose == T) {
        print(paste("rep", i, "completed, final maxZ=", NewMaxZ))
      } 
    }#end of each replicate loop
    return(
      list(
        simmat = simmat,
        alphamat = alphamat,
        indexmat = indexmat,
        phase1Pop = phase1Pop,
        phase2Pop = phase2Pop
      )
    )
  }


# Figure 1: Basic Price Equation and simulations --------------------------

#Panel A
alpha = 4
beta = 1
Max_Z = 1 #arbitrarily set a value for the plot
Npop = 100
Gdist <- rgumbel(Npop, loc=(Max_Z-alpha), scale=beta)
if (new_windows == T){
  quartz()
}
plot(density(Gdist),main="Fig 1A",ylab = expression(Probability~of~matching~italic(z)),xlab=
       "", sub=paste("z"), font.sub=2,bty="n",cex=10, las=1)
abline(v=Max_Z,lty=5,col="red",lwd=2) #Max individual skill
abline(v=Max_Z-alpha,lty=5,lwd=2) #Average learner's attempt

#Panel B
Nbig = 100
Normdist <- rnorm(Nbig,0,1)
big_gen1_df<- as.data.frame(cbind(Normdist,rep(1,length(Normdist))))
big_gen2_df <- as.data.frame(cbind(rgumbel(Nbig, loc=(max(Normdist)-alpha), scale=beta),rep(2,length(Normdist))))
names(big_gen2_df) <- names(big_gen1_df)

big_df<- rbind(big_gen1_df,big_gen2_df) #combine generation 1 and 2
big_df$gen <- as.factor(big_df$V2)
p <- ggplot(big_df, aes(x=gen,y=Normdist))+ 
  geom_violin()+geom_violin()+
  geom_point(aes(x=jitter(V2, factor = .2)),alpha=0.5)+
  geom_point(aes(y=max(big_gen1_df[,1]),x=1),pch=21,bg="red", size = 3)+
  geom_point(aes(y=max(big_gen2_df[,1]),x=2),pch=21,bg="red", size = 3)+
  geom_point(aes(y=max(big_gen1_df[,1]),x=2),pch=1, col="red", size = 3)+
  labs(size=18,x="Generation",y=expression(Skill~level~italic(z)))+ylim(-5,8)+
  theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), axis.line = element_line(colour = "black"),  axis.text=element_text(size=12))


Npop = 20
small_gen1_df <- as.data.frame(cbind(rgumbel(Npop, loc=(max(rnorm(Npop,0,1))-alpha), scale=beta),rep(1,Npop)))
small_gen2_df <- cbind(rgumbel(Npop, loc=(max(small_gen1_df$V1)-alpha), scale=beta),rep(2,Npop))
newsmall <- rbind(small_gen1_df,small_gen2_df)
colnames(newsmall)<-c("Normdist","V2")

newsmall$gen <- as.factor(newsmall$V2)
p2 <- ggplot(newsmall, aes(x=gen,y=Normdist)) + geom_violin()+geom_point(aes(x=jitter(V2,factor = 0.2)),alpha=0.5)+
  geom_point(aes(y=max(newsmall[,1]),x=1),pch=21,bg="red", size = 3)+
  geom_point(aes(y=max(small_gen2_df[,1]),x=2),pch=21,bg="red", size = 3)+
  geom_point(aes(y=max(newsmall[,1]),x=2),pch=1, col="red", size = 3)+
  labs(size=18,x="Generation",y=expression(Skill~level~italic(z)))+ylim(-5,8)+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.line = element_line(colour = "black"),
  axis.text=element_text(size=12)
)
if (new_windows == T){
  quartz()
}
plot_grid(p,p2, labels = "Fig 1B", label_x = 0.05)

#panel C

alphas <- seq(3.1,8,by=0.1)
popsizes <- seq(10,505,by=10)

DZvalmat <- matrix(data=0,nrow = length(alphas),ncol = length(popsizes))

for (i in 1:length(alphas)){
  for (j in 1:length(popsizes)){
    DZvalmat[i,j] <- -alphas[i] + beta*(Euler+log(popsizes[j]))
  }
}

rownames(DZvalmat) <- alphas
colnames(DZvalmat) <- popsizes
abcombs <- expand.grid(rownames(DZvalmat),colnames(DZvalmat))
zvals <- as.list(DZvalmat)
ggvals <- data.frame(abcombs,as.numeric(zvals))
colnames(ggvals) <- c("Alpha","Pop","ΔZ")

panelC<-ggplot(ggvals,aes(x=Pop, y=Alpha,fill=	ΔZ))+geom_tile(col="white")+theme_bw(base_size = 14)+ scale_fill_gradient2(low = "darkred",mid = "white",high = "darkblue")+ggtitle(paste0("Fig 1C: Change per generation in z bar"))+scale_x_discrete(breaks=popsizes[seq(1,length(popsizes),5)])+scale_y_discrete(breaks=alphas[seq(1,length(alphas),10)])+ labs(x= "Population size",y="Variation in \u03b1",cex=5)
if (new_windows == T){
  quartz()
}
panelC

#Panel D

betas <- seq(0.1,5,by=0.1)

for (i in 1:length(betas)){
  for (j in 1:length(popsizes)){
    DZvalmat[i,j] <- -alpha + betas[i]*(Euler+log(popsizes[j]))
  }
}

colnames(DZvalmat) <- popsizes
rownames(DZvalmat) <- betas
abcombs <- expand.grid(rownames(DZvalmat),colnames(DZvalmat))
zvals <- as.list(DZvalmat)
ggvals <- data.frame(abcombs,as.numeric(zvals))
colnames(ggvals) <- c("Beta","Pop","ΔZ")
panelD<-ggplot(ggvals,aes(x=Pop, y=Beta,fill=	ΔZ))+geom_tile(col="white")+theme_bw(base_size = 14)+ scale_fill_gradient2(low = "darkred",mid = "white",high = "darkblue")+ggtitle(paste0("Fig 1D: Change per generation in z bar"))+scale_x_discrete(breaks=popsizes[seq(1,length(popsizes),5)])+scale_y_discrete(breaks=betas[seq(0,length(betas),5)])+ labs(x= "Population size",y="Variation in learning attempts (\u03b2)",cex=5)
if (new_windows == T){
  quartz()
}
panelD


# Figure 2 panel A----------------------------------------------------------------

alpha=5
beta=1
Npop=100
stripereps = 20
paletteLength= stripereps
stripegens = 60
#Set up a custom color palette from red (negative) to darkblue (positive), with white at zero
stripemat <- matrix(0,nrow=stripegens,ncol=stripereps)
Zresults <- rep(0,stripegens)
Max_Z <- 5 #reset the starting Z here for each replicate

#simulate mean Z over time for each replicate
for (i in 1:stripereps){
  NewMaxZ <- Max_Z
  for(j in 1:stripegens){
    Newdist <-  rgumbel(Npop, loc=(NewMaxZ-alpha), scale=beta)  
    NewMaxZ<-max(Newdist)
    meanZ <- mean(Newdist)
    stripemat[j,i] <-meanZ
  }
}


# ed-blue color palette

myColor <- colorRampPalette(c('red','white','darkblue'))(paletteLength)
myBreaks <- c(seq(min(stripemat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(stripemat)/paletteLength, max(stripemat), length.out=floor(paletteLength/2)))

#Order the replicates by where they ended up, from lowest to highest Z value. 
#Find which replicate is the closest positive value to zero -- this is the break index
colrankmat <- stripemat[,order(stripemat[stripegens,], decreasing=F)]
if(max(colrankmat[stripegens,])<0){
  breakind <- stripereps #if they're all negative, the break index is the max value (number of reps)
}else{
  breakind <- which(colrankmat[stripegens,]>0)[1]#otherwise, it's the replicate with lowest positive value 
  
}
breakval <- colrankmat[stripegens,breakind] #the end value of the last generation of this lowest positive value (this will be set to white)
myColor <- colorRampPalette(c('red','white','darkblue'))(paletteLength)
myBreaks <- c(seq(min(colrankmat), breakval, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(colrankmat)/paletteLength, max(colrankmat), length.out=floor(paletteLength/2))) #Make pallette-length string of values to assign colors to, first half (plus one) is from minimum to break value, second half is from (maximum/length) to maximum value
if(max(colrankmat[stripegens,])<0){
  otherind <- 1 #If all the values are negative, "blue" index is the first value
}else{
  otherind<-(stripereps-breakind)-1  #otherwise, blue index is one before the break ind
}

redpal <- colorRampPalette(c("#FF0000", "#FFFCFC"))
bluepal <- colorRampPalette(c("#FCFCFD","#00008B"))
if(max(colrankmat[stripegens,])<0){
  zeropal <- redpal(stripegens) #If all values are negative, just use red pallette for all values
}else{
  zeropal <- c(redpal(breakind),"#FFFFFF",bluepal(otherind)) #The lowest value line is always going to be red, even if it is positive
}




avgstripes <- apply(stripemat,1,mean)

genylim <- c(-50, 50) #adjust if needed
if (new_windows == T){
  quartz()
}
plot(avgstripes,ylim = genylim,pch=21,bg="red",las=1,xlab="generations",ylab="Z (average indicated with open circles)",main=paste("Fig 2A: Alpha=",alpha,", beta=",beta,", N=",Npop,sep=""))
for (i in 1:stripereps){
  linecol <- zeropal[i]
  lines(colrankmat[,i], col= linecol)
  
  lines(lowess(colrankmat[,i],f=.05), col= linecol)
}

points(avgstripes,pch=21,bg="white")
clip(1,100,-100,100)
abline(h=0, lty=2,lwd=2)


# Figure 2 panel B --------------------------------------------------------

beta = 1 
Npop = 100

stripegens=60
syllnum = 5
syll_alphs = c(3,3,3,5,6)
stripemat <- matrix(0,nrow=stripegens,ncol=syllnum)
maxmat <- matrix(0,nrow=stripegens,ncol=syllnum)
reset_Z <- 5 #reset the Z after every replicate
set.seed(234) #to replicate figure in the manuscript
for (i in 1:syllnum){
  curralph <- syll_alphs[i]
  NewMaxZ <- reset_Z
  for (j in 1:stripegens){
    Newdist <-  rgumbel(Npop, loc=(NewMaxZ-curralph), scale=beta)  
    NewMaxZ<-max(Newdist)
    meanZ <- mean(Newdist)
    print(meanZ)
    #Zresults[j] <- NewMaxZ
    stripemat[j,i] <-meanZ
    maxmat[j,i] <- NewMaxZ
  }
  #print(paste("Syll",i,"done"))
}


myColor <- colorRampPalette(c('red','white','gray'))(stripegens)

#Make a palette that's red shades for positive numbers, white for zero, and gray for everything below zero
graycols <- length(which(stripemat[]< 0))
redcols <- length(which(stripemat[]> 0))
redpal <- colorRampPalette(c('white','red'))(redcols)
graypal <- rep("#606060",graycols)
combopal <- append(graypal,redpal)
if (new_windows == T){
  quartz()
}
pheatmap(stripemat,cluster_rows = F, cluster_cols = F,border_color=NA,color = combopal, main = "Fig 2B")





# Figure 3 ----------------------------------------------------------------
# Degree, K, shuffled versus static connections


Npopvals <- c(100,200,300,400)
Kvals <- seq(from=5,to=100,by=10)
reps=3
gens=5
alph = 3
bet = 1

#Panel A. Difference in final Z when each individual samples K individuals in a larger population (N), versus when the total population size is equal to K.
#The difference becomes larger as N gets bigger relative to K, as expected.

finalZ <- rep(0,reps)
Kresults<- matrix(0,nrow = length(Npopvals), ncol = length(Kvals))
for(i in 1:length(Npopvals)){
  Npop = Npopvals[i]
  print(paste("pop",Npop,"started"))
  for(j in 1:length(Kvals)){
    K <- Kvals[j]
    if (K>=Npop){
      finalZ[h] <- NA
    } else {
      for (h in 1:reps){
        fullGraf <- make_full_graph(n=Npop,directed = F,loops = F)
        results <- simulateZ.K(alph=alph,bet=bet,gens=gens,K,fullGraf = fullGraf)
        finalZ[h] <- results$meanZ
      }
    }
    Kresults[i,j] <- mean(finalZ)
  }
  print(paste(Npopvals[i],"Pop done"))
}



Npops <-  seq(from=5,to=100,by=10)
finalZ <- rep(0,reps)
Nresults <- rep(0,length(Npops))


for (h in 1:length(Npops)){
  Popsize <- Npops[h]
  for (i in 1:reps){
    #print(paste("starting pop size",Popsize))
    smallgraf <- make_full_graph(n=Popsize,directed = F,loops = F)
    results <- simulateZ.N(alph=alph,bet=bet,gens=gens,smallgraf=smallgraf,Npop=Popsize)
    finalZ[i] <- results$meanZ
  }
  Nresults[h] <- mean(finalZ)
}

NDegdiff <- sweep(Kresults,2,Nresults) 
if (new_windows == T){
  quartz()
}
pheatmap(NDegdiff[order(nrow(NDegdiff):1),],cluster_rows = F, cluster_cols = F,border_color=NA,color = colorRampPalette(brewer.pal(7,"Blues"))(100), angle_col=0,main=paste("Fig 3A: alpha =",alph))


#Panel B, static versus shuffled
Npop = 100
gens = 10
Dvals <- seq(from=5,to=100,by=10)
Degresults <- rep(0, length(Dvals))

finalZ <- rep(0,reps)
for(j in 1:length(Dvals)) {
  Degval = Dvals[j]
  #print(paste("Starting deg=",Degval))
  Dnet <- sample_degseq(rep(Degval, Npop), method = "vl") #makes network with uniform degree of Degval (K)
  for (h in 1:reps) {
    results <- simulateZ.D(
      alph,
      zStart=1,
      bet,
      gens = gens,
      degnet = Dnet,
      Npop = Npop
    )
    finalZ[h] <- results$meanZ
  }
  #print(paste("i+j=",i+j))
  Degresults[j] <- mean(finalZ)
  #print(paste(Dvals[j], "Deg done"))
}

Kresults_1<- rep(0, length(Kvals))
Npop = 100
for (j in 1:length(Kvals)) {
  K <- Kvals[j]
  #print(paste("Current K", K, "Npop=", Npop))
  
  if (K >= Npop) {
    finalZ[h] <- NA
  } else {
    for (h in 1:reps) {
      fullGraf <- make_full_graph(n = Npop,
                                  directed = F,
                                  loops = F)
      results <-
        simulateZ.K(
          alph = alph,
          bet = bet,
          gens = gens,
          K,
          fullGraf = fullGraf
        )
      finalZ[h] <- results$meanZ
    }
  }
  #print(paste("i+j=",i+j))
  Kresults_1[j] <- mean(finalZ)
  #print(paste("i=",i,"j=",j,"mean last Z=",mean(finalZ)))
 
}


diffmatrix <- data.frame(outer(Degresults,Kresults_1,"-"))
colnames(diffmatrix) <- c("5","15","25","35","45","55","65","75","85","95")

rownames(diffmatrix) <- c("5","15","25","35","45","55","65","75","85","95")
if (new_windows == T){
  quartz()
}
pheatmap(diffmatrix,cluster_rows = F, cluster_cols = F,border_color=NA,color = colorRampPalette(brewer.pal(7,"RdBu"))(100), angle_col=0, main = "Fig3 B")

#Panel C, violin plots 

reps = 5
gens = 50
violinresultsN26 <- data.frame(rep(0,reps),rep("N26"))
colnames(violinresultsN26) = c("meanZ","nettype")

Popsize = 26
for (i in 1:reps){
  smallgraf <- make_full_graph(n=Popsize,directed = F,loops = F)
  results<- simulateZ.N(alph=3,bet=1,gens=gens,smallgraf=smallgraf,Npop = Popsize)
  violinresultsN26[i,1] <- results$meanZ
}

Npop = 100
K=25
violinresults25shuf <- data.frame(rep(0,reps),rep("K 25 Shuffled"))
colnames(violinresults25shuf) = c("meanZ","nettype")

for (h in 1:reps) {
  fullGraf <- make_full_graph(n = Npop,
                              directed = F,
                              loops = F)
  results<-
    simulateZ.K(
      alph = 3,
      bet = 1,
      gens = gens,
      K=K,
      fullGraf = fullGraf
    )
  violinresults25shuf[h,1]<-results$meanZ
}


violinresults25stat <- data.frame(rep(0,reps),rep("K 25 Static"))
colnames(violinresults25stat) = c("meanZ","nettype")

Dnet <- sample_degseq(rep(25, Npop), method = "vl")
for (h in 1:reps) {
  results <- simulateZ.D(
    alph=3,
    zStart=1,
    bet=1,
    gens = gens,
    degnet = Dnet,
    Npop = 100
  )
  #print(paste("results$meanZ",results$meanZ))
  violinresults25stat[h,1] <- results$meanZ
  #print(paste("rep", h, "finished"))
}

plotresults <- rbind(violinresultsN26,violinresults25stat,violinresults25shuf)
p <- ggplot(plotresults,aes(x=nettype,y=meanZ))+geom_violin(draw_quantiles =0.5)+theme_bw() + labs(title = "Fig 3C")
if (new_windows == T){
  quartz()
}
p

# Panel D: see separate network figure script



# Figure 4 ----------------------------------------------------------------

minAlpha= -1 #change to "1" for figure S5.
medAlpha = 5
maxAlpha = 6
simreps = 10

StartZ <- 3
alpha1 <- 5

minThreshZ = 2
maxThreshZ = 10

phase1Pop <- 200
phase1gens <- 500
phase2Pop <- 20 
phase2gens <- 200

BottlePops <-
  c(
    rep(phase1Pop, phase1gens),
    rep(phase2Pop, phase2gens)
  ) #Vector of the population sizes for each generation
simgens <- length(BottlePops)


Newdist <- 0
Gdist <-rgumbel(phase1Pop, loc = (StartZ - alpha1), scale = beta) #Initial probability of matching StartZ
Fig4_results <- AlphaBottleneckSim(StartZ=StartZ, alpha=alpha1,beta=bet, simreps=simreps,alpharange = alpharange, phase1Pop=phase1Pop,phase2Pop=phase2Pop,phase1gens=phase1gens,phase2gens=phase2gens,verbose=F)


simmat <- rbind(rep(StartZ,simgens),Fig4_results$simmat)
avgsims <- apply(simmat,1,mean)
avgalphas <- apply(Fig4_results$alphamat,1,mean)
if (min(simmat) < 0){
  genylim <- c(round(min(simmat))-1, round(max(simmat)+1))
} else {
  genylim <- c(0, round(max(simmat)+1))
}
if (new_windows == T){
  quartz()
}
plot(avgsims,ylim=genylim,pch=".",bg="red",las=1,xlab="generations",ylab=expression(Avg~Population~italic(bar("z"))),main=paste("Fig 4: alpha range=",round(min(Fig4_results$alphamat),2),"-",round(max(Fig4_results$alphamat),2)," (possible:",minAlpha,"-",maxAlpha,") \n beta=",beta," reps=",simreps,", Pops:",phase1Pop,"->",phase2Pop,"at gen",phase1gens))

for (i in 1:simreps){
  
  lines(simmat[,i], col="lightsteelblue")
}

lines(avgsims)
abline(h=0, col="lightgray")
abline(v=phase1gens, col="red")


# Supplemental Figures ----------------------------------------------------

#S1
#Panel A: Varying alphas and K

Npopvals <- c(100)
alphavals <- c(2,3,4,5,6)
Kvals <- seq(from=5,to=100,by=10)
reps=1
gens=5
finalZ <- rep(0,reps)
Kresults<- matrix(0,nrow = length(alphavals), ncol = length(Kvals))
Npop = Npopvals
for(i in 1:length(alphavals)){
  alphaval=alphavals[i]

  for(j in 1:length(Kvals)){

    K <- Kvals[j]
    if (K>=Npop){
      finalZ[h] <- NA
    } else {
      
      for (h in 1:reps){
        fullGraf <- make_full_graph(n=Npop,directed = F,loops = F)
        results <- simulateZ.K(alph=alphaval,bet=1,gens=gens,K,fullGraf = fullGraf)
        finalZ[h] <- results$meanZ
      }
    }
    Kresults[i,j] <- mean(finalZ)
    
  }
  
}
varyingAlpha <- t(Kresults)

max_abs <- max(abs(varyingAlpha))
brk <- do.breaks(c(-max_abs, max_abs), 100)
first_true <- which.max(brk > min(varyingAlpha))
brk <- brk[(first_true -1):length(brk)]
cols <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
cols <- cols[(first_true -1):length(cols)]
if (new_windows == T){
  quartz()
}
pheatmap(varyingAlpha,cluster_rows = F, cluster_cols = F,border_color=NA,color = colorRampPalette(brewer.pal(7,"RdBu"))(100), angle_col=0,main="Fig S1A: varying alphas")
#Panel B: Varying beta
Npopvals <- c(100)
betavals <- c(0.25,0.5,1,1.5,2)
Kvals <- seq(from=5,to=100,by=10)
reps=1
gens=5
finalZ <- rep(0,reps)
Kresults<- matrix(0,nrow = length(betavals), ncol = length(Kvals))
Npop = Npopvals
for(i in 1:length(betavals)){
  betaval=betavals[i]
  
  for(j in 1:length(Kvals)){

    K <- Kvals[j]
    if (K>=Npop){
      finalZ[h] <- NA
    } else {
      
      for (h in 1:reps){
        fullGraf <- make_full_graph(n=Npop,directed = F,loops = F)
        results <- simulateZ.K(alph=3,bet=betaval,gens=gens,K,fullGraf = fullGraf)
        finalZ[h] <- results$meanZ
      }
    }
    Kresults[i,j] <- mean(finalZ)
    
  }
  
}
varyingBeta <- t(Kresults)
max_abs <- max(abs(varyingBeta))
brk <- do.breaks(c(-max_abs, max_abs), 100)
first_true <- which.max(brk > min(varyingBeta))
brk <- brk[(first_true -1):length(brk)]
cols <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
cols <- cols[(first_true -1):length(cols)]
if (new_windows == T){
  quartz()
}
pheatmap(varyingBeta,cluster_rows = F, cluster_cols = F,border_color=NA,color = cols, angle_col=0,main="Fig S1B: varying Betas")


#S3 A. Constant alpha

phase1Pop <- 200
phase1gens <- 500
phase2Pop <- 20
phase2gens <- 500
StartZ <- 1
currAlpha <- 5

BottlePops <-
  c(
    rep(phase1Pop, phase1gens),
    rep(phase2Pop, phase2gens)
  ) #Vector of the population sizes for each generation

simgens <- length(BottlePops)
Zresults <- rep(0, simgens)
simmat <- matrix(0, nrow = simgens, ncol = simreps)

for (i in 1:simreps){
  Newdist <- Gdist
  NewMaxZ <- StartZ
  for(j in 1:simgens){
    currPop <- BottlePops[j]
    
    Newdist <-rgumbel(currPop, loc=(NewMaxZ-currAlpha), scale=beta)  
    NewMaxZ<-max(Newdist) #update max Z
    NewmeanZ <- mean(Newdist)  #update mean Z
    Zresults[j] <- NewMaxZ 
    simmat[j,i] <-NewmeanZ
    
  }
}

simmat <- rbind(rep(StartZ,simgens),simmat)
avgsims <- apply(simmat,1,mean)

#set graph limits
genylim <- c(min(simmat)-1,max(simmat)+1)
if (new_windows == T){
  quartz()
}

plot(avgsims,ylim=genylim,pch=".",bg="red",las=1,xlab="generations",ylab="avg",main=paste("Fig S3A: alpha =",currAlpha,"\n beta=",beta," reps=",simreps,", Pops:",phase1Pop,"->",phase2Pop,"at gen",phase1gens))

for (i in 1:simreps){
  
  lines(simmat[,i], col="lightsteelblue")
}

lines(avgsims)
abline(h=0, col="lightgray")
abline(v=phase1gens, col="red")

#S3 B. Alpha in three "bands"
floor_a <- -1 #if absolute value of floor is <= middle a, trait maintains?
middle_a <- 5
ceiling_a <- 6


alpharange <- c(rep(floor_a,33),rep(middle_a,34),rep(ceiling_a,33)) 
bottleresults <- AlphaBottleneckSim(StartZ=StartZ, alpha=alpha1,beta=bet, simreps=simreps,alpharange = alpharange, phase1Pop=phase1Pop,phase2Pop=phase2Pop,phase1gens=phase1gens,phase2gens=phase2gens,verbose=F)

avgsims <- apply(bottleresults$simmat,1,mean)
maxsims <- apply(bottleresults$simmat,1,max)
avgalphas <- apply(bottleresults$alphamat,1,mean)
avgindex <- apply(bottleresults$indexmat,1,mean)

#set graph limits
if (min(bottleresults$simmat) < 0){
  genylim <- c(round(min(bottleresults$simmat))-1, round(max(bottleresults$simmat)+1))
} else {
  genylim <- c(0, round(max(bottleresults$simmat)+1))
}

if (new_windows == T){
  quartz()
}

plot(avgsims,ylim=genylim,pch=".",bg="red",las=1,xlab="generations",ylab=expression(Avg~Population~italic(bar("z"))),main=paste("Fig S3B: alpha range=",round(min(bottleresults$alphamat),2),"-",round(max(bottleresults$alphamat),2)," (possible:",minAlpha,"-",maxAlpha,") \n beta=",beta," reps=",simreps,", Pops:",bottleresults$phase1Pop,"->",bottleresults$phase2Pop,"at gen",phase1gens))
for (i in 1:simreps){
  
  lines(bottleresults$simmat[,i], col="lightsteelblue")
}

lines(avgsims)
abline(v = phase1gens)
abline(h=0, col="lightgray")


#S4. 
PermZ <- 5
beta <- 1
simreps <- 10
phase1gens <- 500
phase2gens <- 500
simgens <- phase1gens+phase2gens
phase1Pop<-200
phase2Pop <- 20
minAlpha <- -1
maxAlpha <- 6
alpharange <- seq(minAlpha,maxAlpha,length.out = 100)
alpha1 <- alpharange[length(alpharange)/2]

StartZ <- 1

bottleresults <- AlphaBottleneckSim(StartZ=StartZ, alpha=alpha,beta=beta, simreps=simreps,alpharange = alpharange, phase1Pop=phase1Pop,phase2Pop=phase2Pop,phase1gens=phase1gens,phase2gens=phase2gens,verbose=F)

avgsims <- apply(bottleresults$simmat,1,mean)
maxsims <- apply(bottleresults$simmat,1,max)
avgalphas <- apply(bottleresults$alphamat,1,mean)
avgindex <- apply(bottleresults$indexmat,1,mean)

#set graph limits
if (min(bottleresults$simmat) < 0){
  genylim <- c(round(min(bottleresults$simmat))-1, round(max(bottleresults$simmat)+1))
} else {
  genylim <- c(0, round(max(bottleresults$simmat)+1))
}

if (new_windows == T){
  quartz()
}
plot(avgsims,ylim=genylim,pch=".",bg="red",las=1,xlab="generations",ylab=expression(Avg~Population~italic(bar("z"))),main=paste("Figure S4: alpha range=",round(min(bottleresults$alphamat),2),"-",round(max(bottleresults$alphamat),2)," (possible:",minAlpha,"-",maxAlpha,") \n beta=",beta," reps=",simreps,", Pops:",bottleresults$phase1Pop,"->",bottleresults$phase2Pop,"at gen",phase1gens))
for (i in 1:simreps){
  
  lines(bottleresults$simmat[,i], col="lightsteelblue")
}
lines(avgsims)
abline(v = phase1gens)
abline(h=0, col="lightgray")

