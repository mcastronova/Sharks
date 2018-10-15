rm(list=ls())                  # Remove all objects from memory

################## Function Winter Survival ##################
#X is phenotype B
# eventually fecundity based survival
WINTER_SURVIVAL <- function(X, Bs = -.5, As = 4.5)
{
  Survival 					     <- As+Bs*X		# Survival
  #print(Survival)
  return(Survival>0)
} # End of winter survival function

################## Function Relative Fitness ##################
# X is phenotype B
# Function to calculate new mean fitness value and size-based fecundity value
# initial size of reproduction should match intercept
REL_FITNESS <- function(X, Bf = 3, As = 5, Af = 2)
{
  #random sample for survival
  # minimum size to reproduce needs to be added - first few iterations should be growing
  Fecundity 					   <- Af+Bf*X 	# Fecundity
  Fecundity[Fecundity<0] <- 0					# Check on sign
  X.Fitness  				     <- Fecundity # Fitness
  mu.fit							 <- mean(X.Fitness)
  return(X.Fitness/mu.fit)
}	# End of relative fitness function

################## Function Simulation ##################
# N is number of individuals per generations
# MaxYear is number of years
# h2 is heritability (matrix)
# mu (2 number vector)
# P is the mutation rate and should be between .01 and .2
# Ag is intercept of growth curve
# Bg is slope of growth curve
SIMULATION <- function(N = 100, MaxYear = 10, h2 = matrix(c(0.5,0.0, 0, .3), nrow = 2, ncol = 2, byrow = T), mu = c(10, 0.5), P = matrix(c(20, 0, 0, 0.1), byrow = T, nrow = 2, ncol = 2), Ag = 5, Bg = -0.5, mut=0.1)
{
  ################# Main program #################
  set.seed(100)                  # Initialize random number generator
  G			<- diag(sqrt(diag(h2)*diag(P)))      # Calculate Additive genetic variance
  #browser()
  diag(h2) <- 1
  G <- G %*% h2 %*% G
  E <- P-G   #environmental variances and covariances
  ROUNDING <- 1e-8    #to make sample work
  
  #set up array
  Individuals <- array(0,dim = c(N,14,MaxYear)) #All individuals saved in matrix
  colnames(Individuals) <- c("ID", "BirthPhen", "GrowPhen", "BirthGen", "GrowGen", "BirthEnv", "GrowEnv", "Sex", "Fitness", "BirthMut", "GrowMut", "Growth", "Alive", "Age")
  Individuals[,"ID",1] <- c(1:(N)) #add unique IDs
  Individuals[,"Sex",1] <- sample(c(0,1), N, replace = T) #set all sexes
  
  ###### First iteration - setting up first parents ######
  # Age
  Individuals[1:N, "Age",1] <- 10 #sets all indivuduals to 10 as a starting population ****CHANGE THIS IN THE FUTURE***
  
  # Generate Genetic and environmental values using normal distribution
  Individuals[,c("BirthGen", "GrowGen"),1] 	  	<- grfx(N, G, warn = F)   # Genetic values for intercept and slope
  
  # Environmental values
  #use Nadiv grfx() for the random environmental value generation
  Individuals[,c("BirthEnv","GrowEnv"),1] <- grfx(N, E, warn = F)  # Environmental values for intercept and slope
  
  # Phenotype values
  # Phenotypic values (mean(average of parents) + genetic intercept + environmental intercept + genetic slope + environmental slope)*age
  #population mean birth size and population mean growth rate
  Individuals[,"BirthPhen",1]
  Individuals[,"GrowPhen",1] 	 <- (Individuals[,"BirthPhen",1] + Individuals[,"BirthGen",1] + Individuals[,"GrowGen",1] + Individuals[,"BirthEnv",1] + Individuals[,"GrowEnv",1])*Individuals[,"Age",1]
  
  # Winter Survival
  Individuals[,"Alive",1] <- WINTER_SURVIVAL(Individuals[,"GrowPhen",1]) 
  
  # Fitness calculation
  Individuals[Individuals[,"Sex",1] == FALSE,"Fitness",1] <- REL_FITNESS(Individuals[Individuals[,"Sex",1] == FALSE,"GrowPhen",1]) 
  Individuals[Individuals[,"Sex",1] == TRUE,"Fitness",1] <- REL_FITNESS(Individuals[Individuals[,"Sex",1] == TRUE,"GrowPhen",1])
  
  # Mutation values
  # We used reflected gamma distribution from Keightley and Hill (1987), shape paramter specified by them
  # rate value gives mutational variance of .0001*Ve
  # choosing number of mutations - pois distribution with mean of .1 haploid genome; if assume 1000 loci, pois mean .005, they have about 5*`0^(-6) mutation rate per locus
  # if choose between .005 and .1`assumes we have between 1000 and 20000 loci
  # Number of mutations using a Poisson distribution
  # scale parameter .01 and .2 gives us the mean mutations per diploid genome assuming between 1000 and 20000 loci per site mutation rate 5*10^(-6)
  # Bersabe and caballero et al. (2016)
  lambda <- mut*N #mean number of mutations in population
  N.mutations    <- rpois(1,lambda)
  rows <- sample(1:N, N.mutations)
  Individuals[rows,"BirthMut",1] <- Individuals[rows,"BirthMut",1] + sample(c(-1,1), 1) *rgamma(1, shape = 0.5, rate = 10)
  rows <- sample(1:N, N.mutations)
  Individuals[rows,"GrowMut",1] <- Individuals[rows,"GrowMut",1] + sample(c(-1,1), 1) *rgamma(1, shape = 0.5, rate = 10)
  
  #copy to next array before adding mutation
  #num_alive = nrow(Individuals[Individuals[,"Alive",1] == TRUE,,1])
  #if (num_alive == 0){
  #  stop("ALL DIED IN FIRST YEAR!!")
  #}
  #Individuals[1:num_alive,,2] <- Individuals[Individuals[,"Alive",1] == TRUE,,1]
  

  ###### Rest of years from parent values ######
  for (Iyear in 2:MaxYear)          # Iterate over years
  {
    #Copy over alive individuals from previous year
    num_alive = nrow(Individuals[Individuals[,"Alive",Iyear-1] == TRUE,,Iyear-1])
    if (num_alive == 0) {
      print(cat("ALL DIED IN ", Iyear, " YEAR!!"))
      return(Individuals[,,Iyear-1])
    } # End of if
    # all live 
    else if (num_alive == N) {
      # Copy over alive from past generation
      Individuals[1:num_alive,,Iyear] <- Individuals[Individuals[,"Alive",Iyear-1] == TRUE,,Iyear-1]
      
      # Environmental values
      Individuals[,"BrithEnv",Iyear] 	  	<- Individuals[,"Genetic",Iyear] + grfx(N, E, warn = F)  # Environmental values
      Individuals[,"GrowEnv",Iyear] 	  	<- Individuals[,"Genetic",Iyear] + grfx(N, E, warn = F)   # Environmental values 
      
      # Phenotypic values
      #two column matrix of parents average genetic values for a and b, 1/2 both additive variance
      Individuals[,"BirthPhen",Iyear] <- Individuals[,"BirthGen",Iyear] + Individuals[,"BirthEnv",Iyear]
      Individuals[,"GrowPhen",Iyear] <- Individuals[,"GrowPhen",Iyear] <- (Individuals[,"BirthPhen",Iyear] + Individuals[,"BirthGen",Iyear] + Individuals[,"GrowGen",Iyear] + Individuals[,"BirthEnv",Iyear] + Individuals[,"GrowEnv",Iyear])*Individuals[,"Age",Iyear]
      
      # Winter Survival
      Individuals[,"Alive",Iyear] <- WINTER_SURVIVAL(Individuals[,"Phenotype",Iyear]) 
      
      # Fitness calculation
      Individuals[which(Individuals[,"Sex",Iyear] == FALSE & Individuals[,"Alive",Iyear] == TRUE),"Fitness",Iyear] <- REL_FITNESS(Individuals[which(Individuals[,"Sex",Iyear] == FALSE & Individuals[,"Alive",Iyear] == TRUE),"GrowPhen",Iyear])
      Individuals[which(Individuals[,"Sex",Iyear] == TRUE & Individuals[,"Alive",Iyear] == TRUE),"Fitness",Iyear] <- REL_FITNESS(Individuals[which(Individuals[,"Sex",Iyear] == TRUE & Individuals[,"Alive",Iyear] == TRUE),"GrowPhen",Iyear])
    } # End of else if
    else {
      #copy over alive from past generation
      Individuals[1:num_alive,,Iyear] <- Individuals[Individuals[,"Alive",Iyear-1] == TRUE,,Iyear-1]
      #rbind(Individuals[Individuals[,"Alive",Iyear-1] == TRUE,,Iyear-1], matrix(0, nrow = N-num_alive, ncol = 9))

      next_num = N*(Iyear-1)+1
      Individuals[(num_alive+1):N,"ID",Iyear] <- c(next_num:(next_num+(N-(num_alive+1)))) #add unique IDs
      
      #Genotypic value
      print(Iyear)
      
      Iyear_moms = Individuals[which(Individuals[,"Sex",Iyear-1] == TRUE & Individuals[,"Alive",Iyear-1] == TRUE),"GrowGen",Iyear-1] + Individuals[which(Individuals[,"Sex",Iyear-1] == TRUE & Individuals[,"Alive",Iyear-1] == TRUE),"GrowMut",Iyear-1]
      Iyear_dads = Individuals[which(Individuals[,"Sex",Iyear-1] == FALSE & Individuals[,"Alive",Iyear-1] == TRUE),"GrowGen",Iyear-1] + Individuals[which(Individuals[,"Sex",Iyear-1] == FALSE & Individuals[,"Alive",Iyear-1] == TRUE),"GrowMut",Iyear-1]
      # offspring value is average of parents plus random deviation; average of parents+0.5*mutation+random deviation
      #pick parents based on fitness, and grab row index from that year (use sample.int function to sample individuals row in previous year) for parents, mother index vector, father index vector
      Individuals[(num_alive+1):N,"BirthGen",Iyear] <- (sample(Iyear_moms, N-num_alive, replace = TRUE, Individuals[which(Individuals[,"Sex",Iyear-1] == TRUE & Individuals[,"Alive",Iyear-1] == TRUE),"Fitness",Iyear-1]+ROUNDING) 
                                                       + sample(Iyear_dads, N-num_alive, replace = TRUE, Individuals[which(Individuals[,"Sex",Iyear-1] == FALSE & Individuals[,"Alive",Iyear-1] == TRUE),"Fitness",Iyear-1]+ROUNDING))/2 + rnorm(N-num_alive,0,sqrt(mut/2))
     
      # Mutation values
      # Number of mutations using a Poisson distribution
      # P is mutation rate
      lambda         <- 2*mut*(N-num_alive) #mean number of mutations in population
      N.mutations    <- rpois(1,lambda)
      rows <- sample(N, N.mutations)
      Individuals[rows,"BirthMut",Iyear] <- Individuals[rows,"BirthMut",Iyear] + sample(c(-1,1), 1) *rgamma(1, shape = 0.5, rate = 10)
      rows <- sample(1:N, N.mutations)
      Individuals[rows,"GrowMut",Iyear] <- Individuals[rows,"GrowMut",Iyear] + sample(c(-1,1), 1) *rgamma(1, shape = 0.5, rate = 10)
      
      # Sex of Offspring
      Individuals[(num_alive+1):N,"Sex",Iyear] <- sample(c(0,1), N-num_alive, replace = T) #set all sexes

      # Environmental values
      Individuals[,c("BirthEnv", "GrowEnv"),Iyear] 	  	<- Individuals[,"BirthGen",Iyear] + grfx(N, E, warn = F)   # Environmental values
      
      # Phenotypic values
      Individuals[,"GrowPhen",Iyear] <- (Individuals[,"BirthPhen",Iyear] + Individuals[,"BirthGen",Iyear] + Individuals[,"GrowGen",Iyear] + Individuals[,"BirthEnv",Iyear] + Individuals[,"GrowEnv",Iyear])*Individuals[,"Age",Iyear]
      
      
      # Growth
      #Individuals[,"Growth", Iyear] <- Individuals[,"Growth", Iyear] + GROWTH(Individuals[,"agi",Iyear], Individuals[,"bgi",Iyear], Individuals[,"Genetic",Iyear]) #growth on genetic value and adds to old growth value
      #Individuals[,"Phenotype",Iyear] <- Individuals[,"Phenotype",Iyear] + Individuals[,"Growth", Iyear]
      
      # Winter Survival
      Individuals[,"Alive",Iyear] <- WINTER_SURVIVAL(Individuals[,"GrowPhen",Iyear]) 
      
      #browser()
      # Fitness calculation
      Individuals[which(Individuals[,"Sex",Iyear] == FALSE & Individuals[,"Alive",Iyear] == TRUE),"Fitness",Iyear] <- REL_FITNESS(Individuals[which(Individuals[,"Sex",Iyear] == FALSE & Individuals[,"Alive",Iyear] == TRUE),"GrowPhen",Iyear])
      Individuals[which(Individuals[,"Sex",Iyear] == TRUE & Individuals[,"Alive",Iyear] == TRUE),"Fitness",Iyear] <- REL_FITNESS(Individuals[which(Individuals[,"Sex",Iyear] == TRUE & Individuals[,"Alive",Iyear] == TRUE),"GrowPhen",Iyear])
    } # End of else
    
    #add one to ages
    Individuals[,"Age",Iyear+1] <- Individuals[,"Age",Iyear]+1
    
  } # End of Iyear loop2
  
  values <-t(apply(Individuals, MARGIN = c(2,3), FUN = mean))
  #plot(values[,1])
  #par(mfrow=c(1,2))
  #plot(apply(Individuals[,"Phenotype",,drop = TRUE], MARGIN = 1, FUN = var))
  #plot(apply(Individuals[,"Genetic",,drop = TRUE], MARGIN = 1, FUN = var))
  #return(mean(Individuals[,"Phenotype",Iyear]))
  #*****plot birth and growth of individuals at different intervals*******
  #*****histogram of ages******
  #*****plot phenotypes*******
  #scale by mean to compare generations ((sqrt(variance)/mean)*100)
  print(Individuals[,,Iyear])
}

(Individuals <- SIMULATION())
