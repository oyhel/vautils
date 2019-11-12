# Outrages hack to avoid global variable bullshit, but keep it in for now
globalVariables(c("hg18genelist","hg19genelist","hg38genelist", "rsid", "%>%",
                  "chromosome","position","GENE","START","STOP","geneSTART",
                  "geneSTOP", "snp", "chr", "bp", "flanking"))

#' Annotates nearest gene for list of SNPs
#'
#' Finds the nearest gene for a list of markers. The input data expects a data
#' frame with headers: rsid, chromosome, position. The list of
#' annotations are obtained from PLINK documentation
#' https://www.cog-genomics.org/plink/1.9/resources/
#'
#' @param data      dataframe with results
#' @param flanking  flanking distance (in kB) on each side (ie. 20kb = spans 40kb)
#'                  defaults to 100kb (spans 200kb)
#' @param build     string specifying build (hg19 (default), hg18, hg38)
#' @param collapse  if TRUE (default) collapses nearest genes into single commasep.
#'                  line. If FALSE returns all annotated genes on separate lines with
#'                  a 'distance' column showing how many bp upstream (+) or downstream(-)
#'                  the SNP is located
#' @param snp       snp column name in input data
#' @param bp        position/bp column name in input data
#' @param chr       chromosome column name in input data
#' @return data frame with nearest gene specified per marker
#' @export

find_nearest_gene <-function(data, flanking=100, build='hg19', collapse=TRUE, snp='rsid', chr='chromosome', bp='position'){

  data <- data # not sure this is needed

  if(build == 'hg18'){
    genelist <- hg18genelist
  }

  if(build == 'hg19'){
    genelist <- hg19genelist
  }

  if(build == 'hg38'){
    genelist <- hg38genelist
  }

  if(flanking != 100){
    flanking = as.numeric(flanking)
  }

  if(snp != 'rsid'){
    data <- data %>% rename(rsid = snp)
  }

  if(chr != 'chromosome'){
    data <- data %>% rename(chromosome = chr)
  }

  if(bp != 'position'){
    data <- data %>% rename(position = bp)
  }

  data<-sqldf::sqldf(sprintf("select A.*,B.* from
              data A left join genelist B
              ON (A.chromosome == B.CHR and
              A.position >= B.START - %1$s and
              A.position <= B.STOP + %1$s)", flanking*1000)
  )
  if (collapse==TRUE){
    data %>%
      dplyr::group_by(rsid, chromosome, position) %>%
      dplyr::summarise(GENES = paste(GENE, collapse=',')) %>%
      data.frame
  } else {
    data <- data %>%
      dplyr::rename(geneSTART = START, geneSTOP = STOP) %>%
      dplyr::select(rsid, chromosome, position, geneSTART,geneSTOP,GENE)

    data$distance <- apply(data,1, FUN=function(x){
      ifelse(
        !is.na(x['GENE']) & x['position'] < x['geneSTART'], -(as.numeric(x['geneSTART']) - as.numeric(x['position'])),
        ifelse(!is.na(x['GENE']) & x['position'] > x['geneSTOP'], as.numeric(x['position']) - as.numeric(x['geneSTOP']),
               ifelse(!is.na(x['GENE']) & (x['position'] > x['geneSTART']) & (x['position'] < x['geneSTOP']),'intergenic',NA))
      )
    })
    data
  }
}



