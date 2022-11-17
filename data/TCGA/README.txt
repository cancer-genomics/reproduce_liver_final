Downloading TCGA cn data via rtcga

library(RTCGA)
releaseDate="2016-01-28"
downloadTCGA(cancerTypes = "LIHC", destDir = ".", date = releaseDate, dataSet = "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3")
