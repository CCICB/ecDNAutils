% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsing_linx_data.R
\name{vis_segments_replace_centromere_and_telomere_characters_with_position}
\alias{vis_segments_replace_centromere_and_telomere_characters_with_position}
\title{Fix vis_segment centromere/telomer}
\usage{
vis_segments_replace_centromere_and_telomere_characters_with_position(
  vis_segments_df
)
}
\arguments{
\item{vis_segments_df}{Result of \strong{read_vis_segments} (data.frame)}
}
\value{
Data.frame containing updated PosStart and PosEnd columns
}
\description{
Replaces vis_segment start/end positions that were 'C' (centromere) or 'T' (telomere) with absolute base positions
}
\examples{
linx_dir = system.file(package="ecDNAutils", "CHP212")
vis_segments_df <- read_vis_segments(linx_dir, "CHP212")
vis_segments_replace_centromere_and_telomere_characters_with_position(vis_segments_df)
}
