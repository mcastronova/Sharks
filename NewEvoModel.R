rm(list=ls())                  # Remove all objects from memory
#profvis({

################## Function Winter Survival ##################
WINTER_SURVIVAL <- function(X, Bs = -.5, As = 4.5) #X is phenotype
{
  Survival 					     <- As+Bs*X		# Survival
  print(Survival)
  return(Survival>0)
}

################## Function Relative Fitness ##################
#keep relative fitness (size based fecundity)
REL_FITNESS <- function(X, Bf = 3, As = 5, Af = 2) # Function to calculate new mean value with vector and survival and fecundity values
{
  #random sample for survival
  Fecundity 					   <- Af+Bf*X 	# Fecundity
  Fecundity[Fecundity<0] <- 0					# Check on sign
  X.Fitness  				     <- Fecundity # Fitness
  mu.fit							 <- mean(X.Fitness)
  #browser()
  return(X.Fitness/mu.fit)
}	# End of relative fitness function

################## Function Growth ##################
GROWTH <- function(X, Ag = 5, Bg = -0.5)
{
  X <- Ag+Bg*X
  #possibly add to phenotype
}

################## Function Simulation ##################
# N is number of individuals per generations
# MaxYear is number of years
# h2 is heritability
# Vp
# mu
# P is the mutation rate and should be between .01 and .2
SIMULATION <- function(N = 100, MaxYear = 1000, h2 = 0.5, Vp = 1, mu = 10, P = 0.1)
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
    Individuals[X,"Mutation",Iyear] <<- Individuals[X,"Mutation",Iyear] + sample(c(-1,1), 1) *rgamma(1, shape = 0.5, rate = 10)
    #rnorm(1, mean = 0, sd = 0.05)
  } # End function
  
  ################# Main program #################
  set.seed(100)                  # Initialize random number generator
  Va			<- Vp*h2               # Calculate Additive genetic variance
  Ve			<- Vp-Va               # Calculate Environmental variance
  ROUNDING <- 1e-8    #to make sample work
  
  #set up array
  Individuals <- array(0,dim = c(N,9,MaxYear)) #All individuals saved in matrix
  colnames(Individuals) <- c("ID", "Phenotype", "Genetic", "Environment", "Sex", "Fitness", "Mutation", "Growth", "Alive")
  Individuals[,"ID",1] <- c(1:(N)) #add unique IDs
  Individuals[,"Sex",1] <- sample(c(0,1), N, replace = T) #set all sexes
  
  ###### First iteration - setting up first parents ######
  # Generate Genetic and environmental values using normal distribution
  Individuals[,"Genetic",1] 	  	<- rnorm(N, mean=mu, sd=sqrt(Va))   # Genetic values
  Individuals[,"Environment",1] 	  	<- rnorm(N, mean=0,  sd=sqrt(Ve))   # Environmental values 
  Individuals[,"Phenotype",1] 	  	<- Individuals[,"Genetic",1] + Individuals[,"Environment",1]  # Phenotypic values
  
  Individuals[,"Alive",1] <- WINTER_SURVIVAL(Individuals[,"Phenotype",1]) 
  # Fitness calculation
  Individuals[Individuals[,"Sex",1] == FALSE,"Fitness",1] <- REL_FITNESS(Individuals[Individuals[,"Sex",1] == FALSE,"Phenotype",1]) 
  Individuals[Individuals[,"Sex",1] == TRUE,"Fitness",1] <- REL_FITNESS(Individuals[Individuals[,"Sex",1] == TRUE,"Phenotype",1])
  
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
  
  #copy to next array before adding mutation
  num_alive = nrow(Individuals[Individuals[,"Alive",1] == TRUE,,1])
  Individuals[1:num_alive,,2] <- Individuals[Individuals[,"Alive",1] == TRUE,,1]
  
  #adds mutation to genetic value
  #Individuals[,"Genetic",1] <- Individuals[,"Genetic",1] + Individuals[,"Mutation",1] # do this after you copy over so gametes don't affect parent

  ###### Rest of generations from parent values ######
  for (Iyear in 2:MaxYear)          # Iterate over generations
  {
    
    num_alive = nrow(Individuals[Individuals[,"Alive",Iyear-1] == TRUE,,Iyear-1])
    if (num_alive == 0) {
      print("All died")
      return(Individuals[,,Iyear-1])
    }
    # all live 
    else if (num_alive == N) {
      #copy over alive from past generation
      Individuals[1:num_alive,,Iyear] <- Individuals[Individuals[,"Alive",Iyear-1] == TRUE,,Iyear-1]
      
      # Environmental values
      Individuals[,"Environment",Iyear] <- rnorm(N, mean=0,  sd=sqrt(Ve))
      
      # Phenotypic values
      Individuals[,"Phenotype",Iyear] <- Individuals[,"Genetic",Iyear] + Individuals[,"Environment",Iyear] + Individuals[,"Growth",Iyear]
      
      # Growth
      Individuals[,"Growth", Iyear] <- Individuals[,"Growth", Iyear] + GROWTH(Individuals[,"Phenotype",Iyear]) #growth on genetic value and adds to old growth value
      Individuals[,"Phenotype",Iyear] <- Individuals[,"Phenotype",Iyear] + Individuals[,"Growth", Iyear]
      
      # Winter Survival
      Individuals[,"Alive",Iyear] <- WINTER_SURVIVAL(Individuals[,"Phenotype",Iyear]) 
      
      # Fitness calculation
      Individuals[which(Individuals[,"Sex",Iyear] == FALSE & Individuals[,"Alive",Iyear] == TRUE),"Fitness",Iyear] <- REL_FITNESS(Individuals[which(Individuals[,"Sex",Iyear] == FALSE & Individuals[,"Alive",Iyear] == TRUE),"Phenotype",Iyear])
      Individuals[which(Individuals[,"Sex",Iyear] == TRUE & Individuals[,"Alive",Iyear] == TRUE),"Fitness",Iyear] <- REL_FITNESS(Individuals[which(Individuals[,"Sex",Iyear] == TRUE & Individuals[,"Alive",Iyear] == TRUE),"Phenotype",Iyear])
      
    }
    else {
    #copy over alive from past generation
    Individuals[1:num_alive,,Iyear] <- Individuals[Individuals[,"Alive",Iyear-1] == TRUE,,Iyear-1]
    #rbind(Individuals[Individuals[,"Alive",Iyear-1] == TRUE,,Iyear-1], matrix(0, nrow = N-num_alive, ncol = 9))
   
    if(Iyear==2) {
      #browser()
    }
    next_num = N*(Iyear-1)+1
    Individuals[(num_alive+1):N,"ID",Iyear] <- c(next_num:(next_num+(N-(num_alive+1)))) #add unique IDs
    
    #Genotypic value
    print(Iyear)
    
    Iyear_moms = Individuals[which(Individuals[,"Sex",Iyear-1] == TRUE & Individuals[,"Alive",Iyear-1] == TRUE),"Genetic",Iyear-1] + Individuals[which(Individuals[,"Sex",Iyear-1] == TRUE & Individuals[,"Alive",Iyear-1] == TRUE),"Mutation",Iyear-1]
    Iyear_dads = Individuals[which(Individuals[,"Sex",Iyear-1] == FALSE & Individuals[,"Alive",Iyear-1] == TRUE),"Genetic",Iyear-1] + Individuals[which(Individuals[,"Sex",Iyear-1] == FALSE & Individuals[,"Alive",Iyear-1] == TRUE),"Mutation",Iyear-1]
    # offspring value is average of parents plus random deviation; average of parents+0.5*mutation+random deviation
    # add setp to put in mutation so that parents aren't affected by gametes
    #pick parents based on fitness, and grab row index from that year (use sample.int function to sample individuals row in previous year) for parents, mother index vector, father index vector
    Individuals[(num_alive+1):N,"Genetic",Iyear] <- (sample(Iyear_moms, N-num_alive, replace = TRUE, Individuals[which(Individuals[,"Sex",Iyear-1] == TRUE & Individuals[,"Alive",Iyear-1] == TRUE),"Fitness",Iyear-1]+ROUNDING) 
                              + sample(Iyear_dads, N-num_alive, replace = TRUE, Individuals[which(Individuals[,"Sex",Iyear-1] == FALSE & Individuals[,"Alive",Iyear-1] == TRUE),"Fitness",Iyear-1]+ROUNDING))/2 + rnorm(N-num_alive,0,sqrt(Va/2))
    #do mutation here for offspring
    # Mutation values
    # Mean number of mutations in population
    lambda         <- 2*P*(N-num_alive) #mutations for offspring
    # Number of mutations using a Poisson distribution
    N.mutations    <- rpois(1,lambda)
    rows <- sample(N, N.mutations) #instead of N do N-num_alive from parental arrays, offspring getting the mutation
    # Apply function MUTATION to generate mutant loci
    sapply(rows, MUTATION, Iyear = Iyear) #sapply() for list of individuals to get mutational values
    
    # Sex of Offspring
    Individuals[(num_alive+1):N,"Sex",Iyear] <- sample(c(0,1), N-num_alive, replace = T) #set all sexes
    
    # Environmental values
    Individuals[,"Environment",Iyear] <- rnorm(N, mean=0,  sd=sqrt(Ve))
    
    # Phenotypic values
    Individuals[,"Phenotype",Iyear] <- Individuals[,"Genetic",Iyear] + Individuals[,"Environment",Iyear] + Individuals[,"Growth",Iyear]
    
    # Growth
    Individuals[,"Growth", Iyear] <- Individuals[,"Growth", Iyear] + GROWTH(Individuals[,"Phenotype",Iyear]) #growth on genetic value and adds to old growth value
    Individuals[,"Phenotype",Iyear] <- Individuals[,"Phenotype",Iyear] + Individuals[,"Growth", Iyear]
    
    # Winter Survival
    Individuals[,"Alive",Iyear] <- WINTER_SURVIVAL(Individuals[,"Phenotype",Iyear]) 
    
    #browser()
    # Fitness calculation
    Individuals[which(Individuals[,"Sex",Iyear] == FALSE & Individuals[,"Alive",Iyear] == TRUE),"Fitness",Iyear] <- REL_FITNESS(Individuals[which(Individuals[,"Sex",Iyear] == FALSE & Individuals[,"Alive",Iyear] == TRUE),"Phenotype",Iyear])
    Individuals[which(Individuals[,"Sex",Iyear] == TRUE & Individuals[,"Alive",Iyear] == TRUE),"Fitness",Iyear] <- REL_FITNESS(Individuals[which(Individuals[,"Sex",Iyear] == TRUE & Individuals[,"Alive",Iyear] == TRUE),"Phenotype",Iyear])
}
  } # End of Iyear loop2
#})
  
  values <-t(apply(Individuals, MARGIN = c(2,3), FUN = mean))
  #plot(values[,1])
  par(mfrow=c(1,2))
  plot(apply(Individuals[,"Phenotype",,drop = TRUE], MARGIN = 1, FUN = var))
  plot(apply(Individuals[,"Genotype",,drop = TRUE], MARGIN = 1, FUN = var))
  #return(mean(Individuals[,"Phenotype",Iyear]))
  print(Individuals[,,Iyear])
}

#opt <- optim(par = c(.05,3), OPTIM_SIM)

(Individuals <- SIMULATION())

