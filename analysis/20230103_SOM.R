
#Load Libraries
library(tidyverse)
library(rstan)
library(shinystan)
library(bayesplot)
library(broom)
library(vegan)
library(proxy)
library(here)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### Import BRUV data
read.csv(here("data","bruv_species_final_052421.csv")) %>% as_tibble() %>% dplyr::select(-X)-> bruv
### Import Seine data
read.csv(here("data","seine_only_final_data_052421.csv")) %>% as_tibble() %>% dplyr::select(-X)-> seine
### Import eDNA data
e_index_all_long <- readRDS(here("e_index_all_long.RDS")) 

#Load BRUV Data
bruv %>% 
  mutate(., PA= if_else(MaxN>0, 1, 0)) %>% 
  dplyr::select(BRUV_Number, Site,Species, PA) %>% 
  pivot_wider(values_from = PA, names_from = Species) %>% 
  replace(is.na(.), 0) %>% 
  pivot_longer(., cols = `Atherinops_affinis`:`Menticirrhus_undulatus`, names_to = "Species", values_to = "PA") %>% 
  pivot_wider(values_from = PA, names_from = BRUV_Number) %>% 
  mutate(N = `1`+`2`+`3`) %>% 
  dplyr::select(-`1`,-`2`,-`3`) %>% 
  mutate(., K=3 ) %>% 
  mutate(SpeciesIdx = match(Species, unique(Species))) -> bruv_prez_som

#create unique identifier for combinations of site-species; for use in hierarchical modeling
SS_q <- unite(data = bruv_prez_som,
              col = SS_q,
              c("Site", "Species")
) %>% pull(SS_q)
bruv_prez_som$SiteSpecies <- match(SS_q, unique(SS_q)) #index for unique site-species combinations
dataset <-  bruv_prez_som
dataset$N <- as.numeric(dataset$N)
dataset$K <- as.numeric(dataset$K)
#dataset$Date <- as.numeric(dataset$Date)
dataset$SpeciesIdx <- as.numeric(dataset$SpeciesIdx)
#dataset$SiteDateSpecies <- as.numeric(dataset$SiteDateSpecies)
dataset$SiteSpecies <- as.numeric(dataset$SiteSpecies)


##Stan Model

##Stan Model
sink(here("Stan_SOM_hierarchical_bruv.stan"))
cat(
  "data{/////////////////////////////////////////////////////////////////////
  int<lower=1> S;    // number of samples (nrow)
  int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
  int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
  int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
  int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
  int<lower=1> K[S];   // number of replicates per site (ncol)
  int<lower=0> N[S]; // number of detections among these replicates
  int z[S];   // integer flag to help estimate psi parameter
  }
  
  parameters{/////////////////////////////////////////////////////////////////////
  real<lower=0,upper=1> psi[Nloc];  //commonness parameter
  real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
  real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
  }
  
  transformed parameters{/////////////////////////////////////////////////////////////////////
  }
  
  model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
  for (i in 1:S){
  z[i] ~ bernoulli(psi[L[i]]);
  p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
  N[i] ~ binomial(K[i], p[i]);
  }; 
  
  //priors
  psi ~ beta(2,2); 
  p11 ~ beta(2,2); 
  p10 ~ beta(1,10);
  }
  
  generated quantities{
  real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015
  
  for (i in 1:S){
  Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  / ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  + (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
  );
  }
  }
  
  ",
  fill=TRUE)
sink()

#####################
#run Stan model
#note this will take a while the first time you run a particular model, because it needs to compile from C++
#####################      
myHierarchicalModel_bruv <- stan(file = here("Stan_SOM_hierarchical_bruv.stan"), 
                            data = list(
                              S = nrow(dataset),
                              Species = dataset$SpeciesIdx,
                              Nspecies = length(unique(dataset$Species)),
                              L = dataset$SiteSpecies,
                              Nloc = length(unique(dataset$SiteSpecies)),
                              K = dataset$K,
                              N = dataset$N,
                              z = ifelse(dataset$N > 0, 1, 0)
                            ), 
                            chains = 10,   #number of chains
                            iter = 10000   #number of iterations per chain
)

myHierarchicalStanResults_bruv <- tidy(tibble(as.data.frame(myHierarchicalModel_bruv)))
saveRDS(myHierarchicalStanResults_bruv, here("myHierarchicalStanResults_bruv"))

dataset %>% dplyr::select(Species,SpeciesIdx) %>% unique() -> species_matcher_bruv

myHierarchicalStanResults_bruv %>% 
  separate(column, into=c("column","SpeciesIdx"), sep="([\\[\\]])") %>% 
  filter(., column %in% c("p11","p10")) -> species_stan_results_bruv
species_stan_results_bruv %>%  View()
species_stan_results_bruv$SiteSpecies %>%  unique()

species_stan_results_bruv$SpeciesIdx <- as.numeric(species_stan_results_bruv$SpeciesIdx)

species_stan_results_bruv %>% 
  left_join(species_matcher_bruv) %>% 
  pivot_wider(., names_from=column, values_from = c(n, mean, sd,median, trimmed ,   mad   ,  min ,  max, range, skew, kurtosis, se)) %>% 
  group_by(Species) %>% 
  mutate(., Sensitivity = mean_p11/(mean_p11+mean_p10),
         Specificity = (1-mean_p10)/(1-mean_p10+1-mean_p11),
         TPR = mean_p11,
         FPR = mean_p10,
         FNR = 1-mean_p11,
         TNR = 1 - mean_p10) -> bruv_sensitivity
bruv_sensitivity %>%  View()

saveRDS(bruv_sensitivity, file=here("bruv_sensitivity.rds"))

myHierarchicalStanResults_bruv %>% 
  separate(column, into=c("column","SiteSpecies"), sep="([\\[\\]])") %>% 
  filter(., column %in% c("Occupancy_prob")) -> species_stan_results_bruv

species_stan_results_bruv$SiteSpecies <- as.numeric(species_stan_results_bruv$SiteSpecies)

dataset %>% 
  dplyr::select(Species, SiteSpecies) %>%  unique()-> species_matcher_bruv_psi

species_stan_results_bruv %>% 
  left_join(species_matcher_bruv_psi) %>% 
  group_by(Species) %>% 
  dplyr::summarise(., mean_occupancy = mean(mean), sd_occupancy = sd(mean)) -> bruv_occupancy

saveRDS(bruv_occupancy, file=here("bruv_occupancy.rds"))


rm(myHierarchicalModel_bruv)
#p11 is true positive rate TPR
#p10 is false positive rate FPR
#False negative rate os 1-TPR
#True negative rate is 1-FPR
# Sensitivity = TPR/(TPR+FPR)
#Specificity = TNR/(TNR +FNR)


#Load seine Data

seine %>% 
  dplyr::select(Site, Haul_number,Species,Count) %>%
  group_by(Site,Haul_number, Species) %>% 
  dplyr::summarise("Count" = sum(Count)) %>% 
  ungroup() %>% 
  dplyr::select(Site, Haul_number) %>% unique() %>% 
  group_by(Site) %>% 
  count() -> K_holder


seine %>% 
  dplyr::select(Site, Haul_number,Species,Count) %>%
  group_by(Site,Haul_number, Species) %>% 
  dplyr::summarise("Count" = sum(Count)) %>% 
  mutate(., PA= if_else(Count>0, 1, 0)) %>% 
  dplyr::select(Haul_number, Site,Species, PA) %>% 
  pivot_wider(values_from = PA, names_from = Species) %>% 
  replace(is.na(.), 0) %>%
  pivot_longer(., cols = `Amphistichus_argenteus`:`Citharichthys_stigmaeus`, names_to = "Species", values_to = "PA") %>% 
  pivot_wider(values_from = PA, names_from = Haul_number, names_prefix = "Tech") %>% ungroup() %>% 
  mutate(ndetections = rowSums(dplyr::select(., starts_with("Tech")),na.rm = TRUE)) %>% 
  dplyr::select(., !starts_with("Tech"))  %>% 
  left_join(K_holder) %>% 
  mutate(SpeciesIdx = match(Species, unique(Species)),
         N= ndetections,
         K= n) -> seine_prez_som
  
#create unique identifier for combinations of site-species; for use in hierarchical modeling
SS_q <- unite(data = seine_prez_som,
              col = SS_q,
              c("Site", "Species")
) %>% pull(SS_q)
seine_prez_som$SiteSpecies <- match(SS_q, unique(SS_q)) #index for unique site-species combinations
dataset <-  seine_prez_som
dataset$N <- as.numeric(dataset$N)
dataset$K <- as.numeric(dataset$K)
#dataset$Date <- as.numeric(dataset$Date)
dataset$SpeciesIdx <- as.numeric(dataset$SpeciesIdx)
#dataset$SiteDateSpecies <- as.numeric(dataset$SiteDateSpecies)
dataset$SiteSpecies <- as.numeric(dataset$SiteSpecies)


##Stan Model
sink(here("Stan_SOM_hierarchical_seine.stan"))
cat(
  "data{/////////////////////////////////////////////////////////////////////
  int<lower=1> S;    // number of samples (nrow)
  int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
  int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
  int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
  int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
  int<lower=1> K[S];   // number of replicates per site (ncol)
  int<lower=0> N[S]; // number of detections among these replicates
  int z[S];   // integer flag to help estimate psi parameter
  }
  
  parameters{/////////////////////////////////////////////////////////////////////
  real<lower=0,upper=1> psi[Nloc];  //commonness parameter
  real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
  real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
  }
  
  transformed parameters{/////////////////////////////////////////////////////////////////////
  }
  
  model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
  for (i in 1:S){
  z[i] ~ bernoulli(psi[L[i]]);
  p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
  N[i] ~ binomial(K[i], p[i]);
  }; 
  
  //priors
  psi ~ beta(2,2); 
  p11 ~ beta(2,2); 
  p10 ~ beta(1,10);
  }
  
  generated quantities{
  real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015
  
  for (i in 1:S){
  Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  / ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  + (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
  );
  }
  }
  
  ",
  fill=TRUE)
sink()

#####################
#run Stan model
#note this will take a while the first time you run a particular model, because it needs to compile from C++
#####################      
myHierarchicalModel_seine <- stan(file = here("Stan_SOM_hierarchical_seine.stan"), 
                                  data = list(
                                    S = nrow(dataset),
                                    Species = dataset$SpeciesIdx,
                                    Nspecies = length(unique(dataset$Species)),
                                    L = dataset$SiteSpecies,
                                    Nloc = length(unique(dataset$SiteSpecies)),
                                    K = dataset$K,
                                    N = dataset$N,
                                    z = ifelse(dataset$N > 0, 1, 0)
                                  ), 
                                  chains = 10,   #number of chains
                                  iter = 10000   #number of iterations per chain
)

myHierarchicalStanResults_seine <- tidy(tibble(as.data.frame(myHierarchicalModel_seine)))
saveRDS(myHierarchicalStanResults_seine, here("myHierarchicalStanResults_seine"))

dataset %>% dplyr::select(Species,SpeciesIdx) %>% unique() -> species_matcher_seine

myHierarchicalStanResults_seine %>% 
  separate(column, into=c("column","SpeciesIdx"), sep="([\\[\\]])") %>% 
  filter(., column %in% c("p11","p10")) -> species_stan_results_seine


species_stan_results_seine$SpeciesIdx <- as.numeric(species_stan_results_seine$SpeciesIdx)

species_stan_results_seine %>% 
  left_join(species_matcher_seine) %>% 
  pivot_wider(., names_from=column, values_from = c(n, mean, sd,median, trimmed ,   mad   ,  min ,  max, range, skew, kurtosis, se)) %>% 
  group_by(Species) %>% 
  mutate(., Sensitivity = mean_p11/(mean_p11+mean_p10),
         Specificity = (1-mean_p10)/(1-mean_p10+1-mean_p11),
         TPR = mean_p11,
         FPR = mean_p10,
         FNR = 1-mean_p11,
         TNR = 1 - mean_p10) -> seine_sensitivity

saveRDS(seine_sensitivity, file=here("seine_sensitivity.rds"))

myHierarchicalStanResults_seine %>% 
  separate(column, into=c("column","SiteSpecies"), sep="([\\[\\]])") %>% 
  filter(., column %in% c("Occupancy_prob")) -> species_stan_results_seine

species_stan_results_seine$SiteSpecies <- as.numeric(species_stan_results_seine$SiteSpecies)

dataset %>% 
  dplyr::select(Species, SiteSpecies) %>%  unique()-> species_matcher_seine_psi

species_stan_results_seine %>% 
  left_join(species_matcher_seine_psi) %>% 
  group_by(Species) %>% 
  dplyr::summarise(., mean_occupancy = mean(mean), sd_occupancy = sd(mean)) -> seine_occupancy

saveRDS(seine_occupancy, file=here("seine_occupancy.rds"))
rm(myHierarchicalModel_seine)
#p11 is true positive rate TPR
#p10 is false positive rate FPR
#False negative rate os 1-TPR
#True negative rate is 1-FPR
# Sensitivity = TPR/(TPR+FPR)
#Specificity = TNR/(TNR +FNR)


#Load eDNA Data
e_index_all_long %>% 
  mutate(., PA= if_else(e_index>0, 1, 0)) %>% 
  dplyr::select(Species,Site, Bio_rep, Tech_rep,PA) -> edna_long_sens

edna_long_sens %>%  
  pivot_wider(values_from = PA, names_from = Tech_rep, values_fill = list ( PA= 0), names_prefix = "Tech") %>% 
  ungroup() %>% 
  mutate(ndetections = rowSums(dplyr::select(., starts_with("Tech")),na.rm = TRUE)) %>% 
  dplyr::select(., !starts_with("Tech"))  %>% 
  pivot_wider(names_from=Bio_rep, values_from = ndetections, names_prefix = "Bio_rep_") %>% 
  mutate(tot_rep = rowSums(dplyr::select(., starts_with("Bio_rep_")),na.rm = TRUE)) %>% 
  unite(Bio_rep_A, Bio_rep_B, Bio_rep_C, col="pattern_presence") %>%
  unite(Species,Site, col="unique_sample_sum.taxonomy",sep = ":") %>% ungroup() -> Pattern.of.presence_sens


edna_long_sens %>% 
  group_by(Site,Bio_rep) %>% 
  summarise(K=n_distinct(Tech_rep)) %>% 
  ungroup() %>% 
  pivot_wider(names_from=Bio_rep, values_from = K, names_prefix = "Bio_rep_") %>% 
  mutate(K_total = rowSums(dplyr::select(., starts_with("Bio_rep_")),na.rm = TRUE)) %>% 
  unite(Bio_rep_A, Bio_rep_B, Bio_rep_C, col="pattern_tech_reps") %>% 
  mutate(K=3)-> pattern_tech_reps_sens

edna_long_sens %>% 
  unite(Species,Site, col="unique_sample_sum.taxonomy",sep = ":", remove = FALSE) %>% 
  left_join(pattern_tech_reps_sens) %>% 
  left_join(Pattern.of.presence_sens,by ="unique_sample_sum.taxonomy") -> distinct_species_reps_sens

distinct_species_reps_sens %>%  
  pivot_wider(values_from = PA, names_from = Tech_rep, values_fill = list (PA = 0), names_prefix = "Tech") %>% 
  mutate(N = rowSums(dplyr::select(., starts_with("Tech")),na.rm = TRUE)) %>% 
  dplyr::select(., !starts_with("Tech")) -> data_tech_summarized_sens

data_tech_summarized_sens %>% 
  dplyr::select(Species,Site,Bio_rep,N,K) -> data_tech_summarized_sens_ready
unique_data <-data_tech_summarized_sens_ready
#create unique identifier for species; for use in hierarchical modeling
#apparently Ryan's STAN model requires species names to be numbers.
SS_species <- unite(data = unique_data,
                    col = SS_species,
                    c("Species")
) %>% pull(SS_species)
unique_data$Species <- match(SS_species, unique(SS_species)) #index for unique site-species combinations

#create unique identifier for combinations of site-date-species; for use in hierarchical modeling
SDS <- unite(data = unique_data,
             col = SDS,
             c("Site","Bio_rep", "Species")
) %>% pull(SDS)
unique_data$SiteDateSpecies <- match(SDS, unique(SDS)) #index for unique site-date-species combinations

#create unique identifier for combinations of site-species; for use in hierarchical modeling
SS <- unite(data = unique_data,
            col = SS,
            c("Site", "Species")
) %>% pull(SS)
unique_data$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations

unique_data$N <- as.numeric(unique_data$N)


##Stan Model
sink(here("Stan_SOM_hierarchical_edna.stan"))
cat(
  "data{/////////////////////////////////////////////////////////////////////
  int<lower=1> S;    // number of samples (nrow)
  int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
  int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
  int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
  int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
  int<lower=1> K[S];   // number of replicates per site (ncol)
  int<lower=0> N[S]; // number of detections among these replicates
  int z[S];   // integer flag to help estimate psi parameter
  }
  
  parameters{/////////////////////////////////////////////////////////////////////
  real<lower=0,upper=1> psi[Nloc];  //commonness parameter
  real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
  real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
  }
  
  transformed parameters{/////////////////////////////////////////////////////////////////////
  }
  
  model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
  for (i in 1:S){
  z[i] ~ bernoulli(psi[L[i]]);
  p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
  N[i] ~ binomial(K[i], p[i]);
  }; 
  
  //priors
  psi ~ beta(2,2); 
  p11 ~ beta(2,2); 
  p10 ~ beta(1,10);
  }
  
  generated quantities{
  real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015
  
  for (i in 1:S){
  Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  / ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  + (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
  );
  }
  }
  
  ",
  fill=TRUE)
sink()

#####################
#run Stan model
#note this will take a while the first time you run a particular model, because it needs to compile from C++
#####################      
myHierarchicalModel_eDNA <- stan(file = here("Stan_SOM_hierarchical_edna.stan"), 
                                 data = list(
                                   S = nrow(unique_data),
                                   Species = unique_data$Species,
                                   Nspecies = length(unique(unique_data$Species)),
                                   L = unique_data$SiteSpecies,
                                   Nloc = length(unique(unique_data$SiteSpecies)),
                                   K = unique_data$K,
                                   N = unique_data$N,
                                   z = ifelse(unique_data$N > 0, 1, 0)
                                 ), 
                                 chains = 1,   #number of chains
                                 iter = 10000   #number of iterations per chain
)

myHierarchicalStanResults_edna <- tidy(tibble(as.data.frame(myHierarchicalModel_eDNA)))
#launch_shinystan(myHierarchicalModel)
saveRDS(myHierarchicalStanResults_edna, file="myHierarchicalStanResults_edna.rds")

myHierarchicalStanResults_edna <- readRDS(file="myHierarchicalStanResults_edna.rds")

data_tech_summarized_sens_ready$number <- match(data_tech_summarized_sens_ready$Species, unique(data_tech_summarized_sens_ready$Species))

unique_data %>%
  left_join(data_tech_summarized_sens_ready, by = c("Species"="number")) %>% 
  ungroup %>% dplyr::select(SiteSpecies,Species =Species.y) %>% unique() -> species_matcher_edna_psi


unique_data %>%
  left_join(data_tech_summarized_sens_ready, by = c("Species"="number")) %>% 
  ungroup %>% dplyr::select(SpeciesIdx=Species,Species =Species.y) %>% unique() -> species_matcher_edna


myHierarchicalStanResults_edna %>% 
  separate(column, into=c("column","SiteSpecies"), sep="([\\[\\]])") %>% 
  filter(., column %in% c("psi")) -> species_stan_results_eDNA

species_stan_results_eDNA %>%  View()
species_stan_results_eDNA$SiteSpecies <- as.numeric(species_stan_results_eDNA$SiteSpecies)

species_stan_results_eDNA %>% 
  left_join(species_matcher_edna_psi) %>% 
  group_by(Species) %>% 
  dplyr::summarise(., mean_occupancy = mean(mean), sd_occupancy = sd(mean)) -> eDNA_sensitivity_psi

saveRDS(eDNA_sensitivity_psi, file=here("eDNA_sensitivity_psi.rds"))


myHierarchicalStanResults_edna %>% 
  separate(column, into=c("column","SiteSpecies"), sep="([\\[\\]])") %>% 
  filter(., column %in% c("Occupancy_prob")) -> species_stan_results_eDNA

species_stan_results_eDNA %>%  View()
species_stan_results_eDNA$SiteSpecies <- as.numeric(species_stan_results_eDNA$SiteSpecies)

species_stan_results_eDNA %>% 
  left_join(species_matcher_edna_psi) %>% 
  group_by(Species) %>% 
  dplyr::summarise(., mean_occupancy = mean(mean), sd_occupancy = sd(mean)) -> eDNA_occupancy

saveRDS(eDNA_occupancy, file=here("eDNA_occupancy.rds"))

myHierarchicalStanResults_edna %>% 
  separate(column, into=c("column","SpeciesIdx"), sep="([\\[\\]])") %>% 
  filter(., column %in% c("p11","p10")) -> species_stan_results_eDNA

species_stan_results_eDNA %>%  View()
species_stan_results_eDNA$SpeciesIdx <- as.numeric(species_stan_results_eDNA$SpeciesIdx)

species_stan_results_eDNA %>% 
  left_join(species_matcher_edna) %>% 
  pivot_wider(., names_from=column, values_from = c(n, mean, sd,median, trimmed ,   mad   ,  min ,  max, range, skew, kurtosis, se)) %>% 
  group_by(Species) %>% 
  mutate(.,   Sensitivity = mean_p11/(mean_p11+mean_p10),
         Specificity = (1-mean_p10)/(1-mean_p10+1-mean_p11),
         TPR = mean_p11,
         FPR = mean_p10,
         FNR = 1-mean_p11,
         TNR = 1 - mean_p10) -> eDNA_sensitivity

saveRDS(eDNA_sensitivity, file=here("eDNA_sensitivity.rds"))


