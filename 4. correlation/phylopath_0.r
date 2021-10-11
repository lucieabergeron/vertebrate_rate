library(phylopath)
library(ape)
setwd("C:/Users/lu_de/Desktop/work bremen/all_species_data/phologeny/phylopath")

# Tree
tree <- read.tree('phylo_uce.nwk.txt')
my_tree = drop.tip(tree, c("Macaca_mulatta","Elephas_maximus", 
                           "Cercocebus_lunulatus", "Muntiacus_reevesi", 
                           "Oryzias_latipes", "Callithrix_jacchus"))

# Data
my_data <- read.csv2("all_sp_new_new.csv")
my_data <- my_data[1:68,]
rownames(my_data) <- my_data$species_name

# Test
my_data_ne = data.frame(my_data$mu, my_data$generation_time,
                        as.factor(my_data$mating_b), 
                        my_data$mass,
                        my_data$lifespan_wild, 
                        as.factor(my_data$offspring_perg),
                        (my_data$maturation_time_mother + my_data$maturation_time_father)/2)

rownames(my_data_ne) <- my_data$species_name
colnames(my_data_ne) <- c("MU", "GT", "MA", "BM", "LS", "OFF", "MT")

models_ne <- define_model_set(
  one = c(MU ~ GT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  two = c(MU ~ GT + MA, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  three = c(MU ~ GT + MA + OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  four = c(MU ~ GT + MA + OFF + LS, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  five = c(MU ~ GT + MA + OFF + LS + BM, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  six = c(MU ~ GT + MA + OFF + LS + BM + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  seven = c(MU ~ MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  eight = c(MU ~ OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  nine = c(MU ~ LS, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  ten = c(MU ~ BM, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  eleven = c(MU ~ OFF + MA, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  twelve = c(MU ~ OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  thirteen = c(MU ~ OFF + MA + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  fourteen = c(MU ~ OFF + MA + LS, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  fifteen = c(MU ~ OFF + MA + LS + GT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT)
)

models_ne <- define_model_set(
  sixtytwo = c(MU ~ OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  fiftynine = c(MU ~ MA + OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  fiftyfive = c(MU ~ MA + LS + OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  thirtyeight = c(MU ~ GT + OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  fourtyeight = c(MU ~ MA + BM + LS + OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  fourtytwo = c(MU ~ LS + OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  thirtyfive = c(MU ~ GT + MA + OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  four = c(MU ~ BM + LS + OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  thirtynine = c(MU ~ GT + OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  fourtythree = c(MU ~ LS + OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT)
)

plot_model_set(models_ne)
result <- phylo_path(models_ne, data = my_data_ne, tree = my_tree, model = 'lambda')
result
s <- summary(result)
plot(s)
best_model <- best(result)
plot(best_model)
average_model <- average(result)
plot(average_model, algorithm = 'mds', curvature = 0)


# NEW no prior

# Test
my_data_ne = data.frame(my_data$mu, my_data$generation_time,
                        as.factor(my_data$mating_b), 
                        my_data$mass,
                        my_data$lifespan_wild, 
                        as.factor(my_data$offspring_perg),
                        (my_data$maturation_time_mother + my_data$maturation_time_father)/2)

rownames(my_data_ne) <- my_data$species_name
colnames(my_data_ne) <- c("MU", "GT", "MA", "BM", "LS", "OFF", "MT")

models_ne <- define_model_set(
  one = c(MU ~ BM),
  two = c(MU ~ BM + LS),
  three = c(MU ~ BM + LS + MT),
  four = c(MU ~ BM + LS + OFF),
  five = c(MU ~ BM + LS + OFF + MT),
  six = c(MU ~ BM + MT),
  seven = c(MU ~ BM + OFF),
  eight = c(MU ~ BM + OFF + MT),
  nine = c(MU ~ GT),
  ten = c(MU ~ GT + BM),
  eleven = c(MU ~ GT + BM + LS),
  twelve = c(MU ~ GT + BM + LS + MT),
  thirteen = c(MU ~ GT + BM + LS + OFF),
  fourteen = c(MU ~ GT + BM + LS + OFF + MT),
  fifteen = c(MU ~ GT + BM + MT),
  sixteen = c(MU ~ GT + BM + OFF),
  seventeen = c(MU ~ GT + BM + OFF + MT),
  eighteen = c(MU ~ GT + LS),
  nineteen = c(MU ~ GT + LS + MT),
  twenty = c(MU ~ GT + LS + OFF),
  twentyone = c(MU ~ GT + LS + OFF + MT),
  twentytwo = c(MU ~ GT + MA),
  twentythree = c(MU ~ GT + MA + BM),
  twentyfour = c(MU ~ GT + MA + BM + LS),
  twentyfive = c(MU ~ GT + MA + BM + LS + MT),
  twentysix = c(MU ~ GT + MA + BM + LS + OFF),
  twentyseven = c(MU ~ GT + MA + BM + MT),
  twentyeight = c(MU ~ GT + MA + BM + OFF),
  twentynine = c(MU ~ GT + MA + BM + OFF + MT),
  thirty = c(MU ~ GT + MA + LS),
  thirtyone = c(MU ~ GT + MA + LS + MT),
  thirtytwo = c(MU ~ GT + MA + LS + OFF),
  thirtythree = c(MU ~ GT + MA + LS + OFF + MT),
  thirtyfour = c(MU ~ GT + MA + MT),
  thirtyfive = c(MU ~ GT + MA + OFF),
  thirtysix = c(MU ~ GT + MA + OFF + MT),
  thirtyseven = c(MU ~ GT + MT),
  thirtyeight = c(MU ~ GT + OFF),
  thirtynine = c(MU ~ GT + OFF + MT),
  fourty = c(MU ~ LS),
  fourtyone = c(MU ~ LS + MT),
  fourtytwo = c(MU ~ LS + OFF),
  fourtythree = c(MU ~ LS + OFF + MT),
  fourtyfour = c(MU ~ MA),
  fourtyfive = c(MU ~ MA + BM),
  fourtysix = c(MU ~ MA + BM + LS),
  fourtyseven = c(MU ~ MA + BM + LS + MT),
  fourtyeight = c(MU ~ MA + BM + LS + OFF),
  fourtynine = c(MU ~ MA + BM + LS + OFF + MT),
  fifty = c(MU ~ MA + BM + MT),
  fiftyone = c(MU ~ MA + BM + OFF),
  fiftytwo = c(MU ~ MA + BM + OFF + MT),
  fiftythree = c(MU ~ MA + LS),
  fiftyfour = c(MU ~ MA + LS + MT),
  fiftyfive = c(MU ~ MA + LS + OFF),
  fiftysix = c(MU ~ MA + LS + OFF + MT),
  fiftyseven = c(MU ~ MA + MT),
  fiftyeight = c(MU ~ MA + OFF),
  fiftynine = c(MU ~ MA + OFF + MT),
  sixty = c(MU ~ MT),
  sixtyone = c(MU ~ OFF),
  sixtytwo = c(MU ~ OFF + MT)
)

plot_model_set(models_ne)
result <- phylo_path(models_ne, data = my_data_ne, tree = my_tree, model = 'lambda')
result
s <- summary(result)
plot(s)
best_model <- best(result)
plot(best_model)
average_model <- average(result)
plot(average_model, algorithm = 'mds', curvature = 0)

positions <- data.frame(name = c("MA", "GT", "LS", "BM", "MT", "OFF", "MU"),
                        x = c(1, 1.5, 2, 2.5, 3.5, 4, 2),
                        y = c(3, 3, 3, 3, 3, 3, 1) )  
plot(average_model, algorithm = 'mds', curvature = 0,manual_layout = positions)





# Same only for mammals
my_data_mammals = my_data[which(my_data$Group=="Mammal"),]
my_tree_mammals = drop.tip(my_tree, 
                           c("Aptenodytes_forsteri","Ara_glaucogularis",
                             "Bubo_scandiacus","Chauna_torquata","Coturnix_japonica",
                             "Cyanistes_caeruleus","Gallus_gallus","Gyps_fulvus",
                             "Larus_argentatus","Larus_marinus","Pelecanus_crispus",
                             "Phoenicopterus_roseus","Platalea_ajaja","Pygoscelis_adeliae",
                             "Rhea_pennata","Saxicola_maurus","Taeniopygia_guttata",
                             "Turdus_merula","Amphiprion_ocellaris","Betta_splendens",
                             "Cynoglossus_semilaevis","Cyprinus_carpio","Larimichthys_crocea",
                             "Paralichthys_olivaceus","Salmo_salar","Syngnathus_scovelli",
                             "Coleonyx_brevis","Eublepharis_macularius","Pogona_vitticeps",
                             "Sphaerodactylus_macrolepis","Thamnophis_sirtalis","Chrysemys_picta"))


# Test
my_data_ne_mammals = data.frame(my_data_mammals$mu, my_data_mammals$generation_time,
                        as.factor(my_data_mammals$mating_b), 
                        my_data_mammals$mass,
                        my_data_mammals$lifespan_wild, 
                        as.factor(my_data_mammals$offspring_perg),
                        (my_data_mammals$maturation_time_mother + my_data_mammals$maturation_time_father)/2)

rownames(my_data_ne_mammals) <- my_data_mammals$species_name
colnames(my_data_ne_mammals) <- c("MU", "GT", "MA", "BM", "LS", "OFF", "MT")

models_ne_mammals <- define_model_set(
  one = c(MU ~ BM),
  two = c(MU ~ BM + LS),
  three = c(MU ~ BM + LS + MT),
  four = c(MU ~ BM + LS + OFF),
  five = c(MU ~ BM + LS + OFF + MT),
  six = c(MU ~ BM + MT),
  seven = c(MU ~ BM + OFF),
  eight = c(MU ~ BM + OFF + MT),
  nine = c(MU ~ GT),
  ten = c(MU ~ GT + BM),
  eleven = c(MU ~ GT + BM + LS),
  twelve = c(MU ~ GT + BM + LS + MT),
  thirteen = c(MU ~ GT + BM + LS + OFF),
  fourteen = c(MU ~ GT + BM + LS + OFF + MT),
  fifteen = c(MU ~ GT + BM + MT),
  sixteen = c(MU ~ GT + BM + OFF),
  seventeen = c(MU ~ GT + BM + OFF + MT),
  eighteen = c(MU ~ GT + LS),
  nineteen = c(MU ~ GT + LS + MT),
  twenty = c(MU ~ GT + LS + OFF),
  twentyone = c(MU ~ GT + LS + OFF + MT),
  twentytwo = c(MU ~ GT + MA),
  twentythree = c(MU ~ GT + MA + BM),
  twentyfour = c(MU ~ GT + MA + BM + LS),
  twentyfive = c(MU ~ GT + MA + BM + LS + MT),
  twentysix = c(MU ~ GT + MA + BM + LS + OFF),
  twentyseven = c(MU ~ GT + MA + BM + MT),
  twentyeight = c(MU ~ GT + MA + BM + OFF),
  twentynine = c(MU ~ GT + MA + BM + OFF + MT),
  thirty = c(MU ~ GT + MA + LS),
  thirtyone = c(MU ~ GT + MA + LS + MT),
  thirtytwo = c(MU ~ GT + MA + LS + OFF),
  thirtythree = c(MU ~ GT + MA + LS + OFF + MT),
  thirtyfour = c(MU ~ GT + MA + MT),
  thirtyfive = c(MU ~ GT + MA + OFF),
  thirtysix = c(MU ~ GT + MA + OFF + MT),
  thirtyseven = c(MU ~ GT + MT),
  thirtyeight = c(MU ~ GT + OFF),
  thirtynine = c(MU ~ GT + OFF + MT),
  fourty = c(MU ~ LS),
  fourtyone = c(MU ~ LS + MT),
  fourtytwo = c(MU ~ LS + OFF),
  fourtythree = c(MU ~ LS + OFF + MT),
  fourtyfour = c(MU ~ MA),
  fourtyfive = c(MU ~ MA + BM),
  fourtysix = c(MU ~ MA + BM + LS),
  fourtyseven = c(MU ~ MA + BM + LS + MT),
  fourtyeight = c(MU ~ MA + BM + LS + OFF),
  fourtynine = c(MU ~ MA + BM + LS + OFF + MT),
  fifty = c(MU ~ MA + BM + MT),
  fiftyone = c(MU ~ MA + BM + OFF),
  fiftytwo = c(MU ~ MA + BM + OFF + MT),
  fiftythree = c(MU ~ MA + LS),
  fiftyfour = c(MU ~ MA + LS + MT),
  fiftyfive = c(MU ~ MA + LS + OFF),
  fiftysix = c(MU ~ MA + LS + OFF + MT),
  fiftyseven = c(MU ~ MA + MT),
  fiftyeight = c(MU ~ MA + OFF),
  fiftynine = c(MU ~ MA + OFF + MT),
  sixty = c(MU ~ MT),
  sixtyone = c(MU ~ OFF),
  sixtytwo = c(MU ~ OFF + MT)
)

plot_model_set(models_ne_mammals)
result <- phylo_path(models_ne_mammals, data = my_data_ne_mammals, tree = my_tree_mammals, model = 'lambda')
result
s <- summary(result)
plot(s)
best_model <- best(result)
plot(best_model)
average_model <- average(result)
plot(average_model, algorithm = 'mds', curvature = 0)

positions <- data.frame(name = c("MA", "GT", "LS", "BM", "MT", "OFF", "MU"),
                        x = c(1, 1.5, 2, 2.5, 3.5, 4, 2),
                        y = c(3, 3, 3, 3, 3, 3, 1) )  
plot(average_model, algorithm = 'mds', curvature = 0,manual_layout = positions)

best_model <- best(result, boot=500)
