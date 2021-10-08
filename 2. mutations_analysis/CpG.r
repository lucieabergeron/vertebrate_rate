################################################
# CtoA
CtoA_CpG=0
CtoA_no_CpG=0

if(file.size(paste0(path, "pos_CtoA.txt"))==0){
}else{
pos <- read.csv(paste0(path, "pos_CtoA.txt"), sep ="\t", header=FALSE)
data <- read.csv(paste0(path, "CpG_all_CtoA_used.txt"), sep ="\t", header=FALSE)
for (i in 1:nrow(pos)){
    raw_id=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2])
    if(data[raw_id,3]=="G"){
        raw_id_n=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2]-1)
        if(data[raw_id_n,3]=="C"){
            CtoA_CpG=CtoA_CpG+1
        }else{CtoA_no_CpG=CtoA_no_CpG+1}
    }else if (data[raw_id,3]=="C"){
        raw_id_n=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2]+1)
        if(data[raw_id_n,3]=="G"){
            CtoA_CpG=CtoA_CpG+1
        }else{CtoA_no_CpG=CtoA_no_CpG+1}
    }
}}

# CtoG
CtoG_CpG=0
CtoG_no_CpG=0

if(file.size(paste0(path, "pos_CtoG.txt"))==0){
}else{
pos <- read.csv(paste0(path, "pos_CtoG.txt"), sep ="\t", header=FALSE)
data <- read.csv(paste0(path, "CpG_all_CtoG_used.txt"), sep ="\t", header=FALSE)
for (i in 1:nrow(pos)){
    raw_id=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2])
    if(data[raw_id,3]=="G"){
        raw_id_n=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2]-1)
        if(data[raw_id_n,3]=="C"){
            CtoG_CpG=CtoG_CpG+1
        }else{CtoG_no_CpG=CtoG_no_CpG+1}
    }else if (data[raw_id,3]=="C"){
        raw_id_n=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2]+1)
        if(data[raw_id_n,3]=="G"){
            CtoG_CpG=CtoG_CpG+1
        }else{CtoG_no_CpG=CtoG_no_CpG+1}
    }
}}


# CtoT
CtoT_CpG=0
CtoT_no_CpG=0

if(file.size(paste0(path, "pos_CtoT.txt"))==0){
}else{
pos <- read.csv(paste0(path, "pos_CtoT.txt"), sep ="\t", header=FALSE)
data <- read.csv(paste0(path, "CpG_all_CtoT_used.txt"), sep ="\t", header=FALSE)
for (i in 1:nrow(pos)){
    raw_id=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2])
    if(data[raw_id,3]=="G"){
        raw_id_n=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2]-1)
        if(data[raw_id_n,3]=="C"){
            CtoT_CpG=CtoT_CpG+1
        }else{CtoT_no_CpG=CtoT_no_CpG+1}
    }else if (data[raw_id,3]=="C"){
        raw_id_n=which(trimws(as.character(data[,1]))==trimws(as.character(pos[i,1])) & data[,2]==pos[i,2]+1)
        if(data[raw_id_n,3]=="G"){
            CtoT_CpG=CtoT_CpG+1
        }else{CtoT_no_CpG=CtoT_no_CpG+1}
    }
}}

result=matrix(c(CtoA_no_CpG,CtoA_CpG,CtoG_no_CpG,CtoG_CpG,CtoT_no_CpG,CtoT_CpG), nrow=2)
colnames(result)=c("CtoA","CtoG","CtoT")
write.table(result, file=paste0(path, "CpG.tab"), row.names = FALSE)
