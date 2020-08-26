#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  require(corrplot)
  require(argparse)
  require(viridisLite)})


reformat_fp = function(tf, fp, include_patientID, controls, dtype) {
  if (controls == "" |
      (!all(is.element(unique(strsplit(controls, ",")[[1]]),
		      c("PoolTumor", "PoolNormal", 
			"BloodPoolNormal", "NTC"))))) {
    stop(paste0(
	"Unrecognized control sample class. ",
	"Only the following control sample classes are allowed: ",
	paste(c("PoolTumor", "PoolNormal", "NTC"), collapse = ",")))
  }

  names(tf) = gsub("^Sample$", "Sample_ID", names(tf))
  if (any(grepl("^bc", names(fp)))) {
    fp_header = names(fp)
    for (i in 1:nrow(tf)) {
      fp_header = gsub(
        paste0(tf[i,]$Barcode, "_"), 
        paste0(tf[i,]$Sample_ID, "_"),
        fp_header)
    }
    names(fp) = fp_header
  }
  names(fp) = gsub("\\.", "\\-", names(fp))
  non_controls = paste(
    gsub("\\.", "\\-", 
         tf[!tf$Class %in% unique(strsplit(controls, ",")[[1]]),
            ]$Sample_ID), collapse = "|")
  fp = fp[,c("Locus", grep(non_controls, names(fp), value = T))]

  if (!is.null(include_patientID)) {
    check_valid_sample(include_patientID, tf, dtype)
    filter_pattern = paste0(paste(strsplit(include_patientID, ",")[[1]], 
			   collapse = ".*MinorAlleleFreq|"), ".*MinorAlleleFreq")
  } else  {
    filter_pattern = "MinorAlleleFreq"
  }
  fp = fp[,c(
      "Locus", sort(c(grep(filter_pattern, names(fp), value = T))))]
  fp = as.data.frame(lapply(fp, function(y) gsub("^$", "0", y)))
  return(fp)
}


check_valid_sample = function(samples, tf, dtype) {
  for (s in strsplit(samples, ",")[[1]]) {
    if (!s %in% tf$Patient_ID) {
      stop(paste0(s, " is not a valid MRN/patientID for ", dtype, " data"))
    }
  }
}


process_data = function(fp, tf, pID, exclude_controls, dtype) {
  files = lapply(c(tf, fp), function(i){
    read.csv(i, sep="\t", header = T, colClasses = "character")})
  processed_fp = reformat_fp(files[[1]], files[[2]], pID, exclude_controls, dtype)
  return(processed_fp)
}


corr_fp = function(source_fp, source_tf, reference_fp, 
		   reference_tf, source_patientID, reference_patientID, 
		   exclude_controls, outfile) {
  
  source_fp = process_data(source_fp, source_tf, source_patientID,
			   exclude_controls, "source")
  
  if (!(is.null(reference_fp) | is.null(source_fp))) {
    reference_fp = process_data(reference_fp, reference_tf, reference_patientID,
				exclude_controls, "reference")

    common_locus = intersect(source_fp$Locus, reference_fp$Locus)

    fp_matrix = as.data.frame(
      cbind(
        sapply(source_fp[which(source_fp$Locus %in% common_locus),
                         c(grep("Locus", names(source_fp), 
                                value = T, invert = T))], as.numeric, na.rm=T),
        sapply(reference_fp[which(reference_fp$Locus %in% common_locus),
                            c(grep("Locus", names(reference_fp), 
                                   value = T, invert = T))], as.numeric, na.rm=T)
      )
    )
  }
  else {
    common_locus = source_fp$Locus
    fp_matrix = as.data.frame(
      sapply(source_fp[which(source_fp$Locus %in% common_locus),
	     c(grep("Locus", names(source_fp), value = T, invert = T))],
      as.numeric, na.rm=T))
  }

  names(fp_matrix) = gsub("\\.", "-", 
                          gsub("^X","", 
                               gsub("_MinorAlleleFreq", "", names(fp_matrix))))
  fpcor = cor(as.matrix(fp_matrix), method = c("pearson"), use = "complete.obs")
  pdf(outfile, width = 20, height = 20)
  corrplot::corrplot(
    fpcor, tl.cex = 1, col = viridisLite::viridis(
      10, alpha = 0.75, begin = 0, end = 1, direction = -1), 
    tl.col = "black")
  whatever = dev.off()
}

#---------------------------------------------------------------------
args = commandArgs(TRUE)
if (length(args) == 0) {
  message('Run fpcor.R --help for list of input arguments.')
  quit()
}

parser = ArgumentParser(description = 'Calculate fp correlation across samples from two batches')
parser$add_argument('-st', '--source_title_file', required = TRUE,
                    help = 'Title file of the source project')
parser$add_argument('-sp', '--source_FP_file', required = TRUE,
                    help = 'Fingerprint SNP file of the source project')
parser$add_argument('-rt', '--reference_title_file', default = NULL, 
		    help = 'Title file of the reference project')
parser$add_argument('-rp', '--reference_FP_file', default = NULL, 
		    help = 'Fingerprint SNP file of the reference project')
parser$add_argument('-ss', '--source_samples', default = NULL,
                    help = 'comma-separated list of MRNs/PatientIDs that should be selected from the source project for analysis')
parser$add_argument('-rs', '--reference_samples', default = NULL,
                    help = 'comma-separated list of MRNs/PatientIDs that should be selected from the reference project for analysis')
parser$add_argument('-o', '--outfile', required = FALSE, default ='FP_heatmap.pdf', help = 'Name of output file. If full path not provided, current working directory will be used as output directory')
parser$add_argument('-ec', '--exclude_controls', default = "PoolTumor,PoolNormal,NTC",
		    help = 'comma-separated list fo control classes to exclude from analysis. Example: \"NTC,PoolNormal\"')

args = parser$parse_args()
if(basename(args$outfile) == args$outfile) {
  args$outfile = paste0(getwd(), "/", args$outfile)
}

corr_fp(args$source_FP_file,
        args$source_title_file,
        args$reference_FP_file,
        args$reference_title_file,
        args$source_samples,
	args$reference_samples,
	args$exclude_controls,
        args$outfile)
