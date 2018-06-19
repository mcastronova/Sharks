rm(list=ls())                  # Remove all objects from memory
#profvis({

################## Function Relative Fitness ##################
REL_FITNESS <- function(X, Bs, Bf, As = 5, Af = 2) # Function to calculate new mean value with vector and survival and fecundity values
{
  #random sample for survival
  Survival 					     <- As-Bs*X		# Survival
  Survival[Survival<0]   <- 0					# Check on sign
  Fecundity 					   <- Af+Bf*X 	# Fecundity
  Fecundity[Fecundity<0] <- 0					# Check on sign
  X.Fitness  				     <- Survival*Fecundity # Fitness
  mu							 <- mean(X.Fitness)
  #browser()
  return(X.Fitness/mu)
}	# End of relative fitness function

OPTIM_SIM <- function(B = c(.05,3)) {
  Bs=B[1]
  if (length(B) < 2) {
    Bf = B[1]
  }
  else {
    Bf=B[2]
  }
  return(-1*tail(SIMULATION(Bs, Bf),1))
}

################## Function Simulation ##################
# N is number of individuals per generations
# MaxYear is number of years
# h2 is heritability
# Vp
# mu
# P is the mutation rate and should be between .01 and .2
SIMULATION <- function(Bs, Bf, N = 1000, MaxYear = 1000, h2 = 0.5, Vp = 1, mu = 1.5, P = 0.1)
{

  ################## Function Mutation ##################
  MUTATION <- function(X, Iyear = 1)
  {
    # Apply mutation by randomly selecting N.mutations
    #write function to pull random number to add - use  for now rnorm() std=.05
    # We used reflected gamma distribution from Keightley and Hill (1987), shape paramter specified by them
    # rate value gives mutational variance of .0001*Ve
    # choosing number of mutations - pois distribution with mean of .1 haploid genome; if assume 1000 loci, pois mean .005, they have about 5*`0^(-6) mutation rate per locus
    # if choose between .005 and .1`assumes we have between 1000 and 20000 loci 
    Individuals[X,6,Iyear] <<- Individuals[X,6,Iyear] + sample(c(-1,1), 1) *rgamma(1, shape = 0.5, rate = 10)
    #rnorm(1, mean = 0, sd = 0.05)
  } # End function
  
################# Main program #################
set.seed(100)                  # Initialize random number generator
Individuals <- array(0,dim = c(N,6,MaxYear)) #All individuals saved in matrix
Individuals[,4,] <- sample(c(0,1), N, replace = T) #set all sexes
Va			<- Vp*h2               # Calculate Additive genetic variance
Ve			<- Vp-Va               # Calculate Environmental variance
ROUNDING <- 1e-8    #to make sample work

###### First iteration - setting up first parents ######
# Generate Genetic and environmental values using normal distribution
GX 	  	<- rnorm(N, mean=mu, sd=sqrt(Va))   # Genetic values
EX 	  	<- rnorm(N, mean=0,  sd=sqrt(Ve))   # Environmental values 
PX 	  	<- GX + EX + Individuals[,6,1]                     # Phenotypic values

# Combine phenotypic and genetic values
Individuals[,1,1] <- PX
Individuals[,2,1] <- GX
Individuals[,3,1] <- EX

# Fitness calculation
Individuals[Individuals[,4,1] == FALSE,5,1] <- REL_FITNESS(Individuals[Individuals[,4,1] == FALSE,1,1], Bs, Bf)
Individuals[Individuals[,4,1] == TRUE,5,1] <- REL_FITNESS(Individuals[Individuals[,4,1] == TRUE,1,1], Bs, Bf)

# Mutation values
# Mean number of mutations in population
lambda         <- P*N
# Number of mutations using a Poisson distribution
# scale parameter .01 and .2 gives us the mean mutations per diploid genome assuming between 1000 and 20000 loci per site mutation rate 5*10^(-6)
# Bersabe and caballero et al. (2016)
N.mutations    <- rpois(1,lambda)
rows <- sample(1:N, N.mutations)
# Apply function MUTATION to generate mutant loci
sapply(rows, MUTATION) #sapply() for list of individuals to get mutational values

Individuals[,2,1] <- Individuals[,2,1] + Individuals[,6,1]

###### Rest of generations from parent values ######
for (Iyear in 2:MaxYear)          # Iterate over generations
{
  #Genotypic value
  Individuals[,2,Iyear] <- (sample(Individuals[Individuals[,4,Iyear-1] == TRUE,2,Iyear-1], N, replace = TRUE, Individuals[Individuals[,4,Iyear-1] == TRUE,5,Iyear-1]+ROUNDING) 
                    + sample(Individuals[Individuals[,4,Iyear-1] == FALSE,2,Iyear-1], N, replace = TRUE, Individuals[Individuals[,4,Iyear-1] == FALSE,5,Iyear-1]+ROUNDING))/2
  + rnorm(N,0,sqrt(Va/2))
  
  # Environmental values
  Individuals[,3,Iyear] <- rnorm(N, mean=0,  sd=sqrt(Ve))
  
  # Phenotypic values
  Individuals[,1,Iyear] <- Individuals[,2,Iyear] + Individuals[,3,Iyear] 
  
  # Fitness calculation
  Individuals[Individuals[,4,Iyear] == FALSE,5,Iyear] <- REL_FITNESS(Individuals[Individuals[,4,Iyear] == FALSE,1,Iyear], Bs, Bf)
  Individuals[Individuals[,4,Iyear] == TRUE,5,Iyear] <- REL_FITNESS(Individuals[Individuals[,4,Iyear] == TRUE,1,Iyear], Bs, Bf)
  
  # Mutation values
  # Mean number of mutations in population
  lambda         <- P*N
  # Number of mutations using a Poisson distribution
  N.mutations    <- rpois(1,lambda)
  rows <- sample(1:N, N.mutations)
  # Apply function MUTATION to generate mutant loci
  sapply(rows, MUTATION, Iyear = Iyear) #sapply() for list of individuals to get mutational values
  
  Individuals[,2,Iyear] <- Individuals[,2,Iyear] + Individuals[,6,Iyear]
} # End of Iyear loop2
#})

values <-t(apply(Individuals, MARGIN = c(2,3), FUN = mean))
#plot(values[,1])
#par(mfrow=c(1,2))
#plot(apply(Individuals[,1,,drop = TRUE], MARGIN = 2, FUN = var), ylim = c(0,.6))
#plot(apply(Individuals[,2,,drop = TRUE], MARGIN = 2, FUN = var), ylim = c(0,.01))
return(mean(Individuals[,1,Iyear]))
}

#opt <- optim(par = c(.05,3), OPTIM_SIM)

(Individuals <- SIMULATION())
