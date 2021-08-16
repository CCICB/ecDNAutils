
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
  assertthat::assert_that(assertthat::is.dir(path = linx_dir))
  assertthat::assert_that(assertthat::is.string(sample))
  file_to_search_for=paste0(linx_dir, "/", sample, ".linx.links.tsv")

  if(verbose) message("Searching for file: ", file_to_search_for)

  assertthat::assert_that(file.exists(file_to_search_for))

  linx_links_df <- read.csv(file_to_search_for, sep="\t", header=TRUE)

  linx_links_df %>%
    dplyr::filter(ecDna == "true") %>%
    dplyr::select(clusterId, chainId) %>%
    dplyr::distinct() %>%
    return()
}
