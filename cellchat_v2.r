.libPaths("/path/to/your/rlib")
library(CellChat,lib="/path/to/your/rlib")
library(Seurat)
S_Combined.harm = readRDS("S_Combined.harmP.rds")
S_Combined.SCT = readRDS("S_Combined.SCTP.rds")
DefaultAssay(S_Combined.harm) = "RNA"
DefaultAssay(S_Combined.SCT) = "SCT"

doCellChat = function(seuratObject,groups,hs_mm,threads){
    cellchat = createCellChat(object=seuratObject,group.by=groups)
    if (tolower(hs_mm) == "mm"){
        CellChatDB = CellChatDB.mouse
        PPI = PPI.mouse
    }else{
        CellChatDB = CellChatDB.human
        PPI = PPI.human
    }
    cellchat@DB=CellChatDB
    cellchat = subsetData(cellchat)
    future::plan("multisession", workers = as.integer(threads))
    cellchat = identifyOverExpressedGenes(cellchat)
    cellchat = identifyOverExpressedInteractions(cellchat)
    cellchat = projectData(cellchat, PPI) #Uses projected data as opposed to assay data, requires raw.use=FALSE in computeCommunProb function below.
    cellchat = computeCommunProb(cellchat, type = "triMean",raw.use=FALSE)
    cellchat = filterCommunication(cellchat, min.cells = 10)
    cellchat = computeCommunProbPathway(cellchat)
    cellchat = aggregateNet(cellchat)
    return(cellchat)
}

doCellChatSummary = function(cellchatobj){
    #write communication summary
    comm_frame = subsetCommunication(cellchatobj)
    write.csv(comm_frame,"Communication_Summary.csv")
    mat = cellchatobj@net$weight
    myl = list()
    groupSize = as.numeric(table(cellchatobj@idents))
    pdf("Interactions_Weights_Minimal.pdf",width=18,height=18)
    par(mfrow = c(5,4), xpd=TRUE)
    for (i in 1:nrow(mat)) {
        #idx=paste0("P",as.character(i))
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()
    return(comm_frame)
}
LRcontribution = function(cellchatobj,pdf_out){
    pathways.show.all <- cellchatobj@netP$pathways
    myl = list()
    for (i in pathways.show.all){
        z <- netAnalysis_contribution(cellchatobj, signaling = i,return.data=TRUE)
        zp = ggplot(data=z$LR.contribution,aes(x=sort(contribution*100),y=str_replace(name,"\\s*-\\s*","+")))+
            geom_bar(stat="identity")+
            xlab("L+R % Contribution")+ylab("L+R Pair")+
            ggtitle(i)+
            theme(panel.grid.minor= element_blank(),plot.title = element_text(hjust = 0.5))
        myl[[i]] = zp
    }
    g=marrangeGrob(grobs=myl,ncol=9,nrow=9)
    ggsave(paste0(pdf_out,".pdf"),plot=g,width=35,height=35,device="pdf")
}

mychat = doCellChat(seuratObject=S_Combined.harm,groups="immgen.pred",hs_mm="mm",threads=2)
saveRDS(mychat,"HarmonyImmgenCellchat.rds")
mychatsum = doCellChatSummary(cellchatobj=mychat)
LRcontribution(mychat,"LR_Contributions")