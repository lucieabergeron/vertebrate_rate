
# List of sample
data_list=vector()
data_list_full=vector()
for(sam in seq(1,nb_sample)){
    data_list=c(data_list, paste0('data_',sam))
    data_list_full=c(data_list_full, paste0('data_full_',sam))
}

# Per individual
table_each=vector()
for(i in data_list){
    U=length(which(get(i)[,6]=="U"))
    P=length(which(get(i)[,6]=="P"))
    M=length(which(get(i)[,6]=="M"))
    ass=P+M
    prop_p=P/ass
    prop_m=M/ass
    line=c(i, U, P, M, ass, prop_p, prop_m)
    table_each=rbind(table_each, line)
}

table_each=as.data.frame(table_each)
colnames(table_each)<- c("sample", "unassigned", "paternal", "maternal", "assigned", "prop_pat", "prop_mat")

######################################################################################################
# Add total for normalization
table_each_full=vector()
P_tot=0
M_tot=0
for(j in data_list_full){
    U=length(which(get(j)[,6]=="U"))
    P=length(which(get(j)[,6]=="P"))
    M=length(which(get(j)[,6]=="M"))
    line=c(i, U, P, M)
    table_each_full=rbind(table_each_full, line)
    P_tot=P_tot+P
    M_tot=M_tot+M
}
table_each_full=as.data.frame(table_each_full)
colnames(table_each_full)<- c("sample", "unassigned", "paternal_tot", "maternal_tot")

# Final table:
poo <- cbind(table_each, table_each_full[,3:4])
write.table(poo, file=paste0(path, "poo_table.txt"))
# Box plot:
# proportion
#poo_p_p <- poo$prop_pat
#poo_p_m <- poo$prop_mat
# values
#poo_p <- poo$paternal
#poo_m <- poo$maternal
# value normalized
poo_n_p <- as.numeric(as.character(poo$paternal))/as.numeric(as.character(poo$paternal_tot))
poo_n_m <- as.numeric(as.character(poo$maternal))/as.numeric(as.character(poo$maternal_tot))
# proportion of normalized value
poo_np_p <- (poo_n_p)/(poo_n_p+poo_n_m)
poo_np_m <- (poo_n_m)/(poo_n_p+poo_n_m)
# infer number assigned of total number of mutation
#poo_i_p <- round((poo$unassigned+poo$assigned)*poo_np_p)
#poo_i_m <- round((poo$unassigned+poo$assigned)*poo_np_m)
##

if(nb_sample>1){
  test_res = t.test(poo_np_p,poo_np_m)
  png(paste0(path, "norm_prop_box.png"), width = 1300, height = 850)
  par(mar=c(9,17,4,2), mgp=c(3, 3, 0))
  boxplot(poo_np_m, poo_np_p, ylim=c(0,1), names=c("",""), ylab="", cex=3, cex.lab=2.5, cex.axis=2.5, col=c("red2","blue"), yaxt="n")
  axis(2, cex.axis=3.5, las=2)
  text(1.5, test_res$estimate[2], round(test_res$estimate[2],3), cex=2.5, font=2)
  text(1.5, test_res$estimate[1], round(test_res$estimate[1],3), cex=2.5, font=2)
  text(0.6, 0.9, paste0("p =", round(test_res$p.value, 3)), cex=2.5, font=2)
  mtext("Proportion of phased", 2, line=12, cex=4)
  mtext(expression(italic(de~novo)~mutations), 2, line=9, cex=4)
  axis(1, at=c(1,2), label=c("Maternal", "Paternal"), cex.axis=4)
  dev.off()
#
}else{
  png(paste0(path, "norm_prop_box.png"), width = 1300, height = 850)
  par(mar=c(9,17,4,2), mgp=c(3, 3, 0))
  boxplot(poo_np_m, poo_np_p, ylim=c(0,1), names=c("",""), ylab="", cex=3, cex.lab=2.5, cex.axis=2.5, col=c("red2","blue"), yaxt="n")
  axis(2, cex.axis=3.5, las=2)
  text(1.5, poo_np_m, poo_np_m, cex=2.5, font=2)
  text(1.5, poo_np_p, poo_np_p, cex=2.5, font=2)
  mtext("Proportion of phased", 2, line=12, cex=4)
  mtext(expression(italic(de~novo)~mutations), 2, line=9, cex=4)
  axis(1, at=c(1,2), label=c("Maternal", "Paternal"), cex.axis=4)
  dev.off()
}

