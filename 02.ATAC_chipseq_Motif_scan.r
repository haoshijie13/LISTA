library(motifmatchr)

bed_to_granges <- function(file){
   df <- read.table(file,
                    header=F,
                    stringsAsFactors=F)
 
   if(length(df) > 6){
      df <- df[,-c(7:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
   header <- c('chr','start','end','id','score','strand')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
   library("GenomicRanges")
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   }
   return(gr)
}

prepareMotifmatchr <- function(genome, motifs){
	res <- list()

	# get the species name and the genome sequence object based on the object
	genomeObj <- genome
	if (!is.element("BSgenome", class(genomeObj))){
		genomeObj <- getGenomeObject(genome)
	}
	spec <- organism(genomeObj)

	# get the motif PWMs
	motifL <- TFBSTools::PWMatrixList()
	if (is.character(motifs)){
		if (is.element("jaspar", motifs)){
			# copied code from chromVAR, but updated the JASPAR version
			opts <- list()
			opts["species"] <- spec
			opts["collection"] <- "CORE"
			# gets the non-redundant set by default
			mlCur <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
			if (!isTRUE(all.equal(TFBSTools::name(mlCur), names(mlCur)))){
				names(mlCur) <- paste(names(mlCur), TFBSTools::name(mlCur), sep = "_")
			} 
			motifL <- c(motifL, TFBSTools::toPWM(mlCur))
		}
		if (is.element("jaspar_vert", motifs)){
			# JASPER for all vertebrate TFBS
			opts <- list()
			opts["tax_group"] <- "vertebrates"
			opts["collection"] <- "CORE"
			# gets the non-redundant set by default
			mlCur <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
			if (!isTRUE(all.equal(TFBSTools::name(mlCur), names(mlCur)))){
				names(mlCur) <- paste(names(mlCur), TFBSTools::name(mlCur), sep = "_")
			} 
			motifL <- c(motifL, TFBSTools::toPWM(mlCur))
		}
		if (is.element("jaspar2016", motifs)){
			motifL <- c(motifL, TFBSTools::toPWM(chromVAR::getJasparMotifs(species=spec)))
		}
		if (is.element("homer", motifs)){
			if (!requireNamespace("chromVARmotifs")) logger.error(c("Could not load dependency: chromVARmotifs"))
			data("homer_pwms")
			motifL <- c(motifL, chromVARmotifs::homer_pwms)
		}
		if (is.element("encode", motifs)){
			if (!requireNamespace("chromVARmotifs")) logger.error(c("Could not load dependency: chromVARmotifs"))
			data("encode_pwms")
			motifL <- c(motifL, chromVARmotifs::encode_pwms)
		}
		if (is.element("cisbp", motifs)){
			if (!requireNamespace("chromVARmotifs")) logger.error(c("Could not load dependency: chromVARmotifs"))
			if (spec == "Mus musculus"){
				data("mouse_pwms_v1")
				motifL <- c(motifL, chromVARmotifs::mouse_pwms_v1)
			} else if (spec == "Homo sapiens"){
				data("human_pwms_v1")
				motifL <- c(motifL, chromVARmotifs::human_pwms_v1)
			} else {
				logger.warning(c("Could not find cisBP annotation for species", spec))
			}
		}
		if (is.element("cisbp_v2", motifs)){
			if (!requireNamespace("chromVARmotifs")) logger.error(c("Could not load dependency: chromVARmotifs"))
			if (spec == "Mus musculus"){
				data("mouse_pwms_v2")
				motifL <- c(motifL, chromVARmotifs::mouse_pwms_v2)
			} else if (spec == "Homo sapiens"){
				data("human_pwms_v2")
				motifL <- c(motifL, chromVARmotifs::human_pwms_v2)
			} else {
				logger.warning(c("Could not find cisBP annotation for species", spec))
			}
		}
		if (length(motifL) < 1) {
			logger.error(c("No motifs were loaded. Unsupported motifs (?) :", motifs))
		}	
	} else if (is.element("PWMatrixList", class(motifs)) || is.element("PFMatrixList", class(motifs))) {
		motifL <- motifs
	} else {
		logger.error(c("unsupported value for motifs:", motifs))
	}	
	res[["genome"]] <- genomeObj
	res[["motifs"]] <- motifL
	return(res)
}


peak = '/mnt/4/liver_project/atac/01.align/male_liver/Tbl1xr1.anno.bed'
peaks <- bed_to_granges(peak)
cisbp_motif <- prepareMotifmatchr("mm10", "cisbp")$motifs  
sel_motif = c('Jun','Egr1','Fos','Cebpd')
mtf_set = c()
for ( f in sel_motif){
    sel_motif_s <- cisbp_motif[grep(f,names(cisbp_motif))] 
    mtf_set = append(mtf_set, sel_motif_s)
}
motif_pos <- matchMotifs(mtf_set, peaks, genome = "mm10", 
                         out = "positions") 
