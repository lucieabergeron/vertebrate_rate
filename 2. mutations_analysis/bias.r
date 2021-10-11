## source activate hmisc
library(hmisc)
data = read.csv("data.txt", sep='\t')

#####################################################
#
subdata_mammals = data[which(data$Group=="Mammal"),]
subdata_mammals_pat = sum(subdata_mammals$paternal)
subdata_mammals_mat = sum(subdata_mammals$maternal)
subdata_mammals_pat/subdata_mammals_mat
binconf(subdata_mammals_pat, (subdata_mammals_pat+subdata_mammals_mat))
# eg. 0.6951364/0.3048636
# Primates
subdata_primates = subdata_mammals[which(subdata_mammals$order=="Primate"),]
subdata_prim_pat = sum(subdata_primates$paternal)
subdata_prim_mat = sum(subdata_primates$maternal)
subdata_prim_pat/subdata_prim_mat
# No callitrix
subdata_primates_noc = subdata_primates[-c(1:6),]
subdata_prim_noc_pat = sum(subdata_primates_noc$paternal)
subdata_prim_noc_mat = sum(subdata_primates_noc$maternal)
subdata_prim_noc_pat/subdata_prim_noc_mat
binconf(subdata_prim_noc_pat, subdata_prim_noc_pat+subdata_prim_noc_mat)
# Callitrix
subdata_primates_c = subdata_primates[c(1:6),]
subdata_prim_c_pat = sum(subdata_primates_c$paternal)
subdata_prim_c_mat = sum(subdata_primates_c$maternal)
subdata_prim_c_pat/subdata_prim_c_mat
##################################################################
## Scandentia
subdata_sc = subdata_mammals[which(subdata_mammals$order=="Scandentia"),]
subdata_sc_pat = sum(subdata_sc$paternal)
subdata_sc_mat = sum(subdata_sc$maternal)
subdata_sc_pat/subdata_sc_mat
## Rodentia
subdata_rod = subdata_mammals[which(subdata_mammals$order=="Rodentia"),]
subdata_rod_pat = sum(subdata_rod$paternal)
subdata_rod_mat = sum(subdata_rod$maternal)
subdata_rod_pat/subdata_rod_mat
#binconf()
##Artiodactyla + Perissodactyla
subdata_art = subdata_mammals[which(subdata_mammals$order=="Artiodactyla" | subdata_mammals$order=="Perissodactyla"),]
subdata_art_pat = sum(subdata_art$paternal)
subdata_art_mat = sum(subdata_art$maternal)
subdata_art_pat/subdata_art_mat
# Cetartiodactyla
subdata_ceta = subdata_art[c(8,10,21,1,2,3,5,6,7,9,11,12),]
subdata_ceta_pat = sum(subdata_ceta$paternal)
subdata_ceta_mat = sum(subdata_ceta$maternal)
subdata_ceta_pat/subdata_ceta_mat
#binconf
# perisso
subdata_peri = subdata_art[c(4,13,14,15,16,17,18,19,20,22),]
subdata_peri_pat = sum(subdata_peri$paternal)
subdata_peri_mat = sum(subdata_peri$maternal)
subdata_peri_pat/subdata_peri_mat
#binconf

# periso
# subdata_art[c(4,20),]
# ruminantia
#subdata_rumi= subdata_art[c(1,2,3,5,6,7,9,11,12),]
# subdata_rumi_pat = sum(subdata_rumi$paternal)
# subdata_rumi_mat = sum(subdata_rumi$maternal)
# subdata_rumi_pat/subdata_rumi_mat



# Carnivora + Chiroptera
##subdata_car = subdata_mammals[which(subdata_mammals$order=="Carnivora"| subdata_mammals$order=="Chiroptera"),]
subdata_car = subdata_mammals[which(subdata_mammals$order=="Carnivora"),]
subdata_car_pat = sum(subdata_car$paternal)
subdata_car_mat = sum(subdata_car$maternal)
subdata_car_pat/subdata_car_mat


car1=subdata_car[c(1:10,12,13,16),]
sum(car1$paternal)/sum(car1$maternal)
car2=subdata_car[c(11,14,15),]
sum(car2$paternal)/sum(car2$maternal)

# Hyracoidea + Dasyuromorphia + Didelphimorphia
subdata_other = subdata_mammals[which(subdata_mammals$order=="Hyracoidea"| subdata_mammals$order=="Dasyuromorphia" | subdata_mammals$order=="Didelphimorphia"),]
subdata_other_pat = sum(subdata_other$paternal)
subdata_other_mat = sum(subdata_other$maternal)
subdata_other_pat/subdata_other_mat


# Birds
subdata_birds = data[which(data$Group=="Bird"),]
subdata_birds_pat = sum(subdata_birds$paternal)
subdata_birds_mat = sum(subdata_birds$maternal)
subdata_birds_pat/subdata_birds_mat
# Galliformes + Anseriformes
subdata_gal = subdata_birds[which(subdata_birds$order=="Galliformes"| subdata_birds$order=="Anseriformes"),]
subdata_gal_pat = sum(subdata_gal$paternal)
subdata_gal_mat = sum(subdata_gal$maternal)
subdata_gal_pat/subdata_gal_mat
# Passeriformes
subdata_pas = subdata_birds[which(subdata_birds$order=="Passeriformes"),]
subdata_pas_pat = sum(subdata_pas$paternal)
subdata_pas_mat = sum(subdata_pas$maternal)
subdata_pas_pat/subdata_pas_mat
# Charadriiformes +Phoenicopteriformes + Sphenisciformes +Pelecaniformes
subdata_p = subdata_birds[which(subdata_birds$order=="Sphenisciformes"| subdata_birds$order=="Pelecaniformes"),]
subdata_p_pat = sum(subdata_p$paternal)
subdata_p_mat = sum(subdata_p$maternal)
subdata_p_pat/subdata_p_mat

subdata_p = subdata_birds[which(subdata_birds$order=="Charadriiformes" | subdata_birds$order=="Phoenicopteriformes"),]
subdata_p_pat = sum(subdata_p$paternal)
subdata_p_mat = sum(subdata_p$maternal)
subdata_p_pat/subdata_p_mat



# Fish
subdata_fish = data[which(data$Group=="Fish"),]
subdata_fish_pat = sum(subdata_fish$paternal)
subdata_fish_mat = sum(subdata_fish$maternal)
subdata_fish_pat/subdata_fish_mat
# Sub
subdata_fish_less = subdata_fish[c(1,2,6,7,8,10,11,12),]
subdata_fish_less_pat = sum(subdata_fish_less$paternal)
subdata_fish_less_mat = sum(subdata_fish_less$maternal)
subdata_fish_less_pat/subdata_fish_less_mat

# Reptiles
subdata_rep = data[which(data$Group=="Reptile" | data$Group=="Turtle"),]
subdata_rep_pat = sum(subdata_rep$paternal)
subdata_rep_mat = sum(subdata_rep$maternal)
subdata_rep_pat/subdata_rep_mat

# Tutles
##subdata_t = data[which(data$Group=="Turtle"),]
##subdata_t_pat = sum(subdata_t$paternal)
##subdata_t_mat = sum(subdata_t$maternal)
##subdata_t_pat/subdata_t_mat
