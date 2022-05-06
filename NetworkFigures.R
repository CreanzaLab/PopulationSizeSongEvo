#Code to accompany Hudson and Creanza 2022, "Modeling how population size drives the evolution of birdsong, a functional cultural trait" https://doi.org/10.1111/evo.14489

## Make  graphs of networks for figure 3 and S2 in Price Eq manuscript

library(RColorBrewer)
#library(gtools)
library(igraph)
library(ggplot2)
library(evd)
library(sfsmisc)


Euler = 0.5772156649
alph =  3# difficulty to learn, against adaptive evolution of trait
bet = 1 #variance in learning attempts, favors adaptive evolution
Npop = 40 #Needs to be even

reps = 1 


g10 <- sample_degseq(rep(10,Npop),method="vl") #Makes a Npop-member network with degree of 10 for all
l=layout_with_fr(g10)
gen0 <- g10
V(g10)$newZ <- 0
V(g10)$oldZ <- 0
g10_gen0<-g10
g10_results <- rep(0,reps)
gens <- 10 


set.seed(2)
for (h in 1:gens){
    for (j in 1:Npop) {
      V(g10)[j]$oldZ <- V(g10)[j]$newZ
      connectednodes <- neighbors(g10, j)
  
      subs_zmax <- max(V(g10)$newZ[connectednodes])
      gumdist <- rgumbel(Npop, loc = subs_zmax - alph, scale = bet) 
      newZ <- sample(gumdist, 1)
      V(g10)[j]$newZ <- newZ
      V(g10)[j]$color <- newZ
    }
    
}


g10_gen10 <-g10

#Continue for 15 more generations
gens = 15

for (h in 1:gens){
  for (j in 1:Npop) {
    V(g10)[j]$oldZ <- V(g10)[j]$newZ
    connectednodes <- neighbors(g10, j)
    subs_zmax <- max(V(g10)$newZ[connectednodes])
    gumdist <- rgumbel(Npop, loc = subs_zmax - alph, scale = bet) 
    newZ <- sample(gumdist, 1)

    V(g10)[j]$newZ <- newZ
    V(g10)[j]$color <- newZ
  }
 
}


g10_gen25 <-g10
allgens <- c( V(g10_gen10)$newZ, V(g10_gen25)$newZ)
 
D10sf <- max(abs(allgens))
 
##Now make degree five networks
g5 <- sample_degseq(rep(5,Npop),method="vl") #Makes a Npop-member network with degree of 5 
V(g5)$newZ <- 0



g5_gen0 <- g5

V(g5)$newZ <- 0

V(g5)$oldZ <- 0
g5_results <- rep(0,reps)
gens <- 10

set.seed(2)
for (h in 1:gens){
    for (j in 1:Npop) {
      V(g5)[j]$oldZ <- V(g5)[j]$newZ
      connectednodes <- neighbors(g5, j)
 
      subs_zmax <- max(V(g5)$newZ[connectednodes])
      gumdist <- rgumbel(Npop, loc = subs_zmax - alph, scale = bet) 
      newZ <- sample(gumdist, 1)
      V(g5)[j]$newZ <- newZ
      V(g5)[j]$color <- newZ
    }
    #print(paste("Average of generation",h,"=",mean(V(g5)$newZ)))
  }
g5_gen10 <- g5

#Continue for 15 more generations
gens = 15
for (h in 1:gens){
  for (j in 1:Npop) {
    V(g5)[j]$oldZ <- V(g5)[j]$newZ
    connectednodes <- neighbors(g5, j)

    subs_zmax <- max(V(g5)$newZ[connectednodes])
    gumdist <- rgumbel(Npop, loc = subs_zmax - alph, scale = bet) 
    newZ <- sample(gumdist, 1)
    V(g5)[j]$newZ <- newZ
    V(g5)[j]$color <- newZ
  }
}

g5_gen25 <- g5
D5allgens <- c( V(g5_gen10)$newZ, V(g5_gen25)$newZ)



D5sf <- max(abs(D5allgens))


my_palette <- colorRampPalette(c("red", "white", "blue"))(n = Npop)

par(mfrow=c(2,3))

plot(g10_gen0, vertex.label=NA, layout=l, vertex.color="white", main=paste0("Generation 0, Degree = 10, a=",alph," b=",bet), sub= paste0("mean Z of population = ",round(mean(V(g10_gen0)$newZ),2)))
node.colors <- (V(g10_gen10)$newZ+D10sf) / (2*D10sf) * Npop
plot(g10_gen10, vertex.label=NA, layout=l, vertex.color=my_palette[node.colors], main=paste0("Generation 10, Degree = 10, a=",alph," b=",bet), sub= paste0("mean Z of population = ",round(mean(V(g10_gen10)$newZ),2)))
node.colors <- (V(g10_gen25)$newZ+D10sf) / (2*D10sf) * Npop
plot(g10_gen25, vertex.label=NA, layout=l, vertex.color=my_palette[node.colors], main=paste0("Generation 25, Degree = 10, a=",alph," b=",bet), sub= paste0("mean Z of population = ",round(mean(V(g10_gen25)$newZ),2)))


plot(g5_gen0, vertex.label=NA, layout=l, vertex.color="white", main=paste0("Generation 0, Degree = 5, a=",alph," b=",bet), sub= paste0("mean Z of population = ",round(mean(V(g5_gen0)$newZ),2)))
node.colors <- (V(g5_gen10)$newZ+D5sf) / (2*D5sf) * Npop
plot(g5_gen10, vertex.label=NA, layout=l, vertex.color=my_palette[node.colors], main=paste0("Generation 10, Degree = 5, a=",alph," b=",bet), sub= paste0("mean Z of population = ",round(mean(V(g5_gen10)$newZ),2)))
node.colors <- (V(g5_gen25)$newZ+D5sf) / (2*D5sf) * Npop
plot(g5_gen25, vertex.label=NA, layout=l, vertex.color=my_palette[node.colors], main=paste0("Generation 25, Degree = 5, a=",alph," b=",bet),sub= paste0("mean Z of population = ",round(mean(V(g5_gen25)$newZ),2)))


