#' run_cytotrace_pipline
#'
#' @param rds      seurat object
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#' @param chunk.size  How many cells to include in each sub-matrix
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'
run_cytotrace_pipline <- function(
    rds,
    assay = NULL,
    slot = "counts",
    ncores = 1,
    chunk.size=1000
){
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }
  matrix <- Seurat::GetAssayData(object = rds,
                                 slot = slot,
                                 assay = assay)
  res <- CytoTRACE(as.matrix(matrix), batch= NULL, enableFast=TRUE, ncores=ncores, subsamplesize=chunk.size )
  rds <- AddMetaData(rds, res$CytoTRACE, col.name='scEntropyscore_cytotrace')
  return(rds)

}

#' run_slice_pipline
#'
#' @param rds      seurat object
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param species species choose,only choose Human, Mouse,Rat
#' @param chunk.size  How many cells to include in each sub-matrix
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'
run_slice_pipline <- function(
    rds,
    assay = NULL,
    slot = "data",
    species = "human",
    chunk.size=1000
){
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }
  matrix <- Seurat::GetAssayData(object = rds,
                                 slot = slot,
                                 assay = assay)

  matrix <- as.data.frame(matrix)
  #  matrix <- matrix[rowSums(matrix) > 0, ]
  #  dim(matrix)

  sc <- construct(exprmatrix=matrix,
                  cellidentity=Idents(rds),
                  projname="SLICE")
  if (!(species %in% c("Human", "Mouse", "Rat"))){
    stop("Species only support Human, Mouse, Rat")
  }
  if (species == "Human"){
    load(system.file("rda", "hs_kappasim.rda", package = "scEntropyR"))
    km <- kh
  }
  if (species == "Mouse"){
    load(system.file("rda", "mm_kappasim.rda", package = "scEntropyR"))
    km <- km
  }
  if (species == "Rat"){
    load(system.file("rda", "rat_kappasim.rda", package = "scEntropyR"))
    km <- km_rat
  }
  sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of mouse genes
                   calculation="bootstrap",               # choose the bootstrap calculation
                   B.num=100,                             # 100 iterations
                   exp.cutoff=1,                          # the threshold for expressed genes
                   B.size=chunk.size,                           # the size of bootstrap sample
                   clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters
                   random.seed=123)
  rds <- AddMetaData(rds, sc@entropies, col.name='scEntropyscore_slice')
  return(rds)
}


#' run_stemid_pipline
#'
#' @param rds      seurat object
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'
run_stemid_pipline <- function(
    rds,
    assay = NULL,
    slot = "counts"
){
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }
  matrix <- Seurat::GetAssayData(object = rds,
                                 slot = slot,
                                 assay = assay)

  matrix <- as.matrix(matrix)
  data <- t(t(matrix)/apply(matrix,2,sum))*median(apply(matrix,2,sum)) + 0.1
  probs <- t(t(data)/apply(data,2,sum))
  score <- data.frame(-apply(probs*log(probs)/log(nrow(data)),2,sum))
  rds <- AddMetaData(rds,score , col.name='scEntropyscore_stemid')
  return(rds)

}

#' run_scent_ccat_pipline
#'
#' @param rds      seurat object
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param species species choose,only choose Human
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'
run_scent_ccat_pipline <- function(
    rds,
    assay = NULL,
    slot = "counts",
    species = "Human"
){
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }
  matrix <- Seurat::GetAssayData(object = rds,
                                 slot = slot,
                                 assay = assay)
  if (!(species %in% c("Human"))){
    stop("Species only support Human")
  }
  if (species == "Human") {
    gene <- mapIds(x=org.Hs.eg.db,keys = rownames(matrix),
                   column='ENTREZID', keytype='SYMBOL', multiVals = "first" )
    gene <- data.frame(unlist(gene))
    gene <- na.omit(gene)
    colnames(gene)<-"id"
    matrix <- matrix[rownames(gene),]
    rownames(matrix)<-gene$id
    score <- SCENT::CompCCAT(exp.m = as.matrix(matrix), ppiA.m = SCENT::net13Jun12.m)
    rds <- AddMetaData(rds, score, col.name='scEntropyscore_scent_ccat')
    return(rds)
  }
}

#' run_scent_sr_pipline
#'
#' @param rds      seurat object
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param species species choose,only choose Human
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'
run_scent_sr_pipline <- function(
    rds,
    assay = NULL,
    slot = "counts",
    species = "Human",
    ncores=1
){
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }
  matrix <- Seurat::GetAssayData(object = rds,
                                 slot = slot,
                                 assay = assay)
  if (!(species %in% c("Human"))){
    stop("Species only support Human")
  }
  if (species == "Human") {
    gene <- mapIds(x=org.Hs.eg.db,keys = rownames(matrix),
                   column='ENTREZID', keytype='SYMBOL', multiVals = "first" )
    gene <- data.frame(unlist(gene))
    gene <- na.omit(gene)
    colnames(gene)<-"id"
    matrix <- matrix[rownames(gene),]
    rownames(matrix)<-gene$id
    integ <- DoIntegPPI(matrix,SCENT::net13Jun12.m)
    sr_res <- CompSRana(integ,local=TRUE,mc.cores =ncores )
    sr_score <- sr_res$SR
    rds <- AddMetaData(rds, sr_score, col.name='scEntropyscore_scent_sr')
    return(rds)
  }
}



#' run_entropy
#'
#' @param rds      seurat object
#' @param method Method selection parameters for calculating cellular entropy，you can choose "cytotrace","slice","stemid","scent_sr" and "scent_ccat"
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#' @param species species choose,When the method is set to "scent_ccat" or "scent_sr", the species parameter can only be selected as "Human". When the method is set to "slice", the species parameter can be selected from Human, Mouse, and Rat. Other methods are not affected by the species parameter.
#' @param chunk.size  How many cells to include in each sub-matrix
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'
run_entropy <- function(
    rds,
    method,
    assay = NULL,
    slot = NULL,
    ncores = 1,
    species = NULL,
    chunk.size=1000
){

  if (method %in% c("cytotrace","slice","stemid" ,"scent_ccat","scent_sr")) {
    print(paste0("------------------run ", method, "------------------"))
    timestart <- Sys.time()
    if (method == "cytotrace"){
      if (is.null(slot)) {
        slot <- "counts"
      }
      rds1 <- run_cytotrace_pipline(rds=rds,assay = assay,slot=slot,ncores = ncores,chunk.size=chunk.size)
    }
    if (method == "slice"){
      if (is.null(slot)) {
        slot <- "data"
      }
      rds1 <- run_slice_pipline(rds=rds,assay = assay,slot=slot,species=species,chunk.size=chunk.size)
    }
    if (method == "stemid"){
      if (is.null(slot)) {
        slot <- "counts"
      }
      rds1 <- run_stemid_pipline(rds=rds,assay = assay,slot=slot)
    }
    if (method == "scent_ccat"){
      if (is.null(slot)) {
        slot <- "counts"
      }
      if (is.null(species)){
        stop("please check species!")
      }
      rds1 <- run_scent_ccat_pipline(rds=rds,assay = assay,slot=slot,species = species)
    }
    if (method == "scent_sr"){
      if (is.null(slot)) {
        slot <- "counts"
      }
      if (is.null(species)){
        stop("please check species!")
      }
      rds1 <- run_scent_sr_pipline(rds=rds,assay = assay,slot=slot,species = species,ncores=ncores)
    }
    timeend <- Sys.time()
    time_use <- round(difftime(timeend, timestart, units = "secs"), 2)
    print(paste0("time use :", time_use, " s"))
    print(paste0(
      "------------------run ",
      method,
      " done!------------------"
    ))
    return(rds1)
  } else {
    stop(paste0("please check method !"))
  }

}
