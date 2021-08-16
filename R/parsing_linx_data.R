
#' identify ecDNA chains
#'
#' Lists clusters/chain combinations that comprise ecDNA molecules in a single sample
#'
#' @param linx_dir path to a directory containing linx files including samplename.linx.links.tsv (string)
#' @param sample name of sample of interest (string)
#' @param verbose verbose (bool)
#'
#' @return data.frame describing SV clusters and chains that comprise ecDNA molecules in a single sample
#' @export
#'
#' @examples
#' linx_dir = system.file(package="utilitybeltlinx", "CHP212")
#' get_ecDNA_chains(linx_dir, "CHP212")
get_ecDNA_chains <- function(linx_dir, sample, verbose = TRUE){

  linx_links_df <- read_links(linx_dir, sample, verbose)

  linx_links_df %>%
    dplyr::filter(ecDna == "true") %>%
    dplyr::select(clusterId, chainId) %>%
    dplyr::distinct() %>%
    return()
}

#' get_ecDNA_chain_coordinates
#'
#' Get coordinates of segments in ecDNA molecules
#'
#' @param linx_dir path to a directory containing linx files
#' @param sample name of sample to examine
#' @param cluster name of clusters of interest.
#' By default will return all clusters/chains in ecDNA molecules.
#' Can be vector of elements, where element 1 of cluster corresponds to element 1 of chainId, element 2 of cluster corresponds to element 2 of chainId ... etc
#' @param chainId name of chain ids of interest.
#' By default will return all clusters/chains in ecDNA molecules.
#' Can be vector of elements, where element 1 of chainId corresponds to element 1 of cluster, element 2 of chainId corresponds to element 2 of cluster ... etc
#' @return coordinates of a single ecDNA molecule (data.frame)
#' @export
#'
#' @examples
#' linx_dir = system.file(package="utilitybeltlinx", "CHP212")
#' get_ecDNA_chain_coordinates(linx_dir, "CHP212", cluster=274, chainId = 8)
get_ecDNA_chain_coordinates <- function(linx_dir, sample, cluster = NA, chainId = NA){
  vis_segments_df=read_vis_segments(linx_dir, sample, verbose = FALSE)
  links_df <- read_links(linx_dir, sample, verbose = FALSE)

  ecdna_cluster_chain_ids <- links_df %>%
    dplyr::filter(ecDna == "true") %>%
    dplyr::mutate(cluster_chain_id = paste(clusterId, chainId)) %>%
    dplyr::pull(cluster_chain_id) %>%
    unique()

  ecDNA_segments <- vis_segments_df %>%
    dplyr::filter(SampleId == sample) %>%
    dplyr::mutate(cluster_chain_id = paste(ClusterId,ChainId)) %>%
    dplyr::filter(cluster_chain_id %in% ecdna_cluster_chain_ids)

  if(!is.na(cluster) & !is.na(chainId)){
    assertthat::assert_that(length(cluster)==length(chainId), msg = "length of cluster and chainID must be the same")
    user_specified_cluster_chain_id = paste(cluster, chainId)
    ecDNA_segments <- ecDNA_segments %>% dplyr::filter(cluster_chain_id %in% user_specified_cluster_chain_id)
  }

  return(ecDNA_segments)
}

#' ecDNA coords to bed
#'
#' Get coordinates of segments in ecDNA molecules and write to bedfile
#'
#' @inheritParams get_ecDNA_chain_coordinates
#' @param bed_outfile path to output bedfile
#'
#' @return bedfile describing ecDNA coordinates
#' @export
#'
get_ecDNA_chain_coordinates_write_bed <- function(linx_dir, sample, cluster, chainId, bed_outfile){
  ecDNA_segments <- get_ecDNA_chain_coordinates(linx_dir, sample, cluster, chainId)
  vis_segments_to_bed_file(ecDNA_segments, bed_outfile)
}


# export ------------------------------------------------------------------

#' vis_segments to bed
#'
#' @param vis_segments_df result of \strong{read_vis_segments} (data.frame)
#'
#' @return bed-like data.frame
#' @export
#'
#' @examples
#' linx_dir = system.file(package="utilitybeltlinx", "CHP212")
#' vis_segments_df <- read_vis_segments(linx_dir, "CHP212")
#' vis_segments_to_bed_df(vis_segments_df)
vis_segments_to_bed_df <- function(vis_segments_df){
  assertthat::assert_that(is.data.frame(vis_segments_df))

  vis_segments_df %>%
    dplyr::mutate(sample_cluster_chain_id = paste(SampleId, ClusterId, ChainId, sep = "__")) %>%
    dplyr::select(Chromosome, PosStart, PosEnd, sample_cluster_chain_id)
}

#' vis_segments to bed
#'
#' Writes a vis_segments dataframe to disk as a bed file
#'
#' @inheritParams vis_segments_to_bed_df
#'
#' @param outfile path to output bedfile
#'
#' @return
#' @export
#'
#' @examples
#' linx_dir = system.file(package="utilitybeltlinx", "CHP212")
#' vis_segments_df <- read_vis_segments(linx_dir, "CHP212")
#' vis_segments_to_bed_file(vis_segments_df)
vis_segments_to_bed_file <- function(vis_segments_df, outfile, verbose=TRUE){
  assertthat::assert_that(assertthat::is.string(outfile))

  if(verbose) message("Writing bed file to: ", outfile)

  vis_segments_to_bed_df(vis_segments_df = vis_segments_df) %>%
    write.table(file = outfile,
                sep = "\t",
                col.names = FALSE,
                quote = FALSE,
                row.names = FALSE
    )
}

# Preprocessing --------------------------------------------------------------

#' Fix vis_segment centromere/telomer
#'
#' Replaces vis_segment start/end positions that were 'C' (centromere) or 'T' (telomere) with absolute base positions
#'
#' @param vis_segments_df result of \strong{read_vis_segments} (data.frame)
#'
#' @return data.frame containing updated PosStart and PosEnd columns
#' @export
#'
#' @examples
#' linx_dir = system.file(package="utilitybeltlinx", "CHP212")
#' vis_segments_df <- read_vis_segments(linx_dir, "CHP212")
#' vis_segments_replace_centromere_and_telomere_characters_with_position(vis_segments_df)
vis_segments_replace_centromere_and_telomere_characters_with_position <- function(vis_segments_df){
  assertthat::assert_that(is.data.frame(vis_segments_df))

  message("Assuming reference genome was hg19/GRCh37")
  message("Replacing start/end positions that were 'C' (centromere) or 'T' (telomere) with absolute base positions")

  chromosomes_hg19_df <- structure(
    list(
      chrom = c(
        "1",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "2",
        "20",
        "21",
        "22",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "X",
        "Y"
      ),
      centromere = c(
        125000000L,
        40200000L,
        53700000L,
        35800000L,
        17900000L,
        17600000L,
        19000000L,
        36600000L,
        24000000L,
        17200000L,
        26500000L,
        93300000L,
        27500000L,
        13200000L,
        14700000L,
        91000000L,
        50400000L,
        48400000L,
        61000000L,
        59900000L,
        45600000L,
        49000000L,
        60600000L,
        12500000L
      ),
      chromLength = c(
        249250621L,
        135534747L,
        135006516L,
        133851895L,
        115169878L,
        107349540L,
        102531392L,
        90354753L,
        81195210L,
        78077248L,
        59128983L,
        243199373L,
        63025520L,
        48129895L,
        51304566L,
        198022430L,
        191154276L,
        180915260L,
        171115067L,
        159138663L,
        146364022L,
        141213431L,
        155270560L,
        59373566L
      )
    ),
    row.names = c(NA,
                  -24L),
    class = c("tbl_df", "tbl", "data.frame")
  )

  segments_in_ecdna_df <- vis_segments_df %>%
    dplyr::mutate(dplyr::across(c(PosStart, PosEnd), .fns = as.character))

  segments_in_ecdna_coordinates_fixed_df <- segments_in_ecdna_df %>%
    dplyr::mutate(
      PosStartOld = PosStart,
      PosStart = sapply(seq_along(PosStart), function(i){
        if(PosStart[i] == "T"){
          return("0")
        }
        else if(PosStart[i] == "C"){
          chromosomes_hg19_df[chromosomes_hg19_df$chrom==Chromosome[i],"centromere"] %>%
            as.character() %>%
            return()
        }
        else
          return(PosStart[i])
      }) %>% as.numeric(),
      PosEndOld = PosEnd,
      PosEnd = sapply(seq_along(PosEnd), function(i){
        if(PosEnd[i] == "T"){
          chromosomes_hg19_df[chromosomes_hg19_df$chrom==Chromosome[i],"chromLength"] %>%
            as.character() %>%
            return()
        }
        else if(PosEnd[i] == "C"){
          chromosomes_hg19_df[chromosomes_hg19_df$chrom==Chromosome[i],"centromere"] %>%
            as.character() %>%
            return()
        }
        else
          return(PosEnd[i])
      }) %>% as.numeric()
    )

  return(segments_in_ecdna_coordinates_fixed_df)
}


# File Reading ------------------------------------------------------------
read_linx_file <- function(linx_dir, sample, suffix, verbose=TRUE){
  assertthat::assert_that(assertthat::is.dir(path = linx_dir))
  assertthat::assert_that(assertthat::is.string(sample))
  file_to_search_for=paste0(linx_dir, "/", sample, suffix)

  if(verbose) message("Searching for file: ", file_to_search_for)
  assertthat::assert_that(file.exists(file_to_search_for))
  if(verbose) message("    > File found")

  read.csv(file_to_search_for, header = TRUE, sep = "\t")
}

read_vis_segments <- function(linx_dir, sample, verbose=TRUE, fix_centromeres_and_telomeres=TRUE){
  vis_segment_df <- read_linx_file(linx_dir = linx_dir, sample = sample, suffix = ".linx.vis_segments.tsv", verbose = verbose)

  if(fix_centromeres_and_telomeres)
    vis_segment_df <- vis_segments_replace_centromere_and_telomere_characters_with_position(vis_segment_df)
}

read_vis_sv_data <- function(linx_dir, sample, verbose=TRUE){
  read_linx_file(linx_dir = linx_dir, sample = sample, suffix = ".linx.vis_sv_data.tsv", verbose = verbose)
}

read_links <- function(linx_dir, sample, verbose=TRUE){
  read_linx_file(linx_dir = linx_dir, sample = sample, suffix = ".linx.links.tsv", verbose = verbose)
}

read_ecdna <- function(linx_dir, sample, verbose=TRUE){
  read_linx_file(linx_dir = linx_dir, sample = sample, suffix = ".linx.ecdna.csv", verbose = verbose)
}
