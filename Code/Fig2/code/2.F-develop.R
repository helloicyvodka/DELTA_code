F.develop <- function(gene.init, gene.network, tmax = 50, Dmax = 10, a = 100, bias = 0){
  gene.OF <- F.geneOF(gene.init)
  cell.tree  <- Node$new(name = "Root", DAge = 0, tAge = 0, Lineage="", Gene.level = gene.init, Gene.OF = gene.OF)
  tAge = 0
  DAge = 0
  while (TRUE){
    if (tAge > tmax || DAge > Dmax) break
    for (cell.leaf in cell.tree$leaves){
      cell.divide = F.divide(cell.leaf$Gene.level)
      if (is.null(cell.divide)){
        tAge  <- cell.leaf$tAge + 1
        if (tAge > tmax || DAge > Dmax) break
        name  <- "LR"
        DAge <- cell.leaf$DAge
        Lineage <- cell.leaf$Lineage
        Gene.level <- F.grow(cell.leaf$Gene.level, gene.network = gene.network, a = a,bias = bias)
        Gene.OF <- F.geneOF(Gene.level)
        #                 cell.leaf$AddChild(name = name, DAge = DAge, tAge = tAge, 
        #                                    Lineage = Lineage, Gene.level = Gene.level, Gene.OF = Gene.OF)
        
        #                 cell.leaf$name <- name
        cell.leaf$DAge <- DAge
        cell.leaf$tAge <- tAge
        cell.leaf$Lineage <- Lineage
        cell.leaf$Gene.level <- Gene.level
        cell.leaf$Gene.OF <- Gene.OF
      }
      else {
        tAge <- cell.leaf$tAge + 1
        DAge <- cell.leaf$DAge + 1
        if (tAge > tmax || DAge > Dmax) break
        Lname <- "L"
        LLineage <- paste0(cell.leaf$Lineage,"0")
        Rname <- "R"
        RLineage <- paste0(cell.leaf$Lineage,"1")
        LRinit <- F.divide(cell.leaf$Gene.level)
        
        LGene.level <- F.grow(LRinit$L, gene.network = gene.network, a = a,bias=bias)
        RGene.level <- F.grow(LRinit$R, gene.network = gene.network, a = a,bias=bias)
        LGene.OF  <- F.geneOF(LGene.level)
        RGene.OF <- F.geneOF(RGene.level)
        cell.leaf$AddChild(name = Lname, DAge = DAge, tAge = tAge, 
                           Lineage = LLineage, Gene.level = LGene.level, Gene.OF = LGene.OF)
        cell.leaf$AddChild(name = Rname, DAge = DAge, tAge = tAge, 
                           Lineage = RLineage, Gene.level = RGene.level, Gene.OF = RGene.OF)
      }
    }
  }
  return (cell.tree)
}