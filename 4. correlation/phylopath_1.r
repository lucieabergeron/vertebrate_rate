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
  one = c(MU ~ GT, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~ GT),
  two = c(MU ~ GT, GT ~ MA, GT ~ BM, GT ~ LS, GT ~ OFF, GT ~ MT),
  three = c(MU ~ OFF, OFF ~ GT, OFF ~ BM, OFF ~ LS, OFF ~ MT, OFF ~ MA),
  four = c(MU ~ BM, BM ~ MA, BM ~ GT, BM ~ LS, BM ~ OFF, BM ~ MT),
  five = c(MU ~ OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~GT),
  six = c(MU ~ OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~GT),
  seven = c(MU ~ OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~GT, LS ~ BM, OFF ~MA),
  eight = c(MU ~ OFF + MT + LS, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~GT, LS ~ BM, OFF ~MA),
  nine = c(MU ~ OFF + MT + BM, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~GT, LS ~ BM, OFF ~MA),
  ten = c(MU ~ OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~GT, LS ~ BM, OFF ~ MA + MT),
  eleven = c(MU ~ OFF + MT + BM, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~GT, LS ~ BM, OFF ~MA),
  twelve = c(MU ~ OFF + MT + LS, MA ~ GT, BM ~ GT, LS ~ GT, OFF ~ GT, MT ~ GT, LS ~ BM, OFF ~ MA),
  thirteen = c(MU ~ GT + OFF, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT),
  fourteen = c(MU ~ MA + OFF + MT, MA ~ GT, BM ~ GT, LS ~ GT, OFF~GT, MT ~GT)
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
                        x = c(3, 3, 4, 4, 2, 2, 1),
                        y = c(1.5, 2, 2.5, 1.5, 3, 1.5, 2) )  
plot(average_model, algorithm = 'mds', curvature = 0,manual_layout = positions)
plot(best_model,manual_layout = positions)

