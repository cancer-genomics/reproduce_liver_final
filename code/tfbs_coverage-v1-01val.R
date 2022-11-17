library(BSgenome.Hsapiens.UCSC.hg19)

test_args <- function() {
  index <<- 1
  fragDir <<- "/dcs04/scharpf/data/annapragada/Liver_Validation/granges"
  outDir <<- "/dcs04/scharpf/data/zfoda/remap202_val"
}

# Getting command line arguments
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])
fragDir <- args[2]
outDir <- args[3]

# Getting the list of files containing fragments
frag.files <- list.files(fragDir)

# Defining the sample name
s <- gsub(".rds$", "", frag.files[index])

# if this sample was already processed, skipping...
if (file.exists(file.path(outDir, "tmp", paste0(s, "_rel_cov.rds")))) {
  stop("This sample was already processed. Quitting...")
}

# Reading in the fragments for sample 'index'
frag.paths <- file.path(fragDir, frag.files)
frag.gr <- readRDS(frag.paths[index])

# Assigning seqlengths to frag.gr
chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
seqinfo(frag.gr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[chrs]

# Computing per base coverage
cov <- coverage(frag.gr)
rm(frag.gr)

# For each set of TFBSs, compute coverage
tfbs.dir <- "/dcl01/scharpf1/data/dbruhm/delfi_followup/tfbs-pipeline/data/individual_tfs"
tfbs.sets <- list.files(tfbs.dir)

# Defining regions for computing a relative coverage metric
ctl_pos <- c(-3000:-2500, 2500:3000)
dip_pos <- -100:100

# rel_cov will store the relative coverage for each TF
rel_cov <- rep(NA, length(tfbs.sets))
for (i in 1:length(tfbs.sets)) {
  print(i)  

  tfbs.gr <- readRDS(file.path(tfbs.dir, tfbs.sets[i]))

  # Computing the mean coverage for each TSS
  # If there are >25,000 TFBSs, splitting to avoid integer overflow
  n.tfbs <- length(tfbs.gr)
  n.chunks <- ceiling(n.tfbs/25000)
  chunk.id <- rep(1:n.chunks, each = 25000)[1:n.tfbs]
  tfbs.gr <- split(tfbs.gr, chunk.id)
  rm(chunk.id)

  # Computing the coverage at each position relative to each TFBS
  cts <- vector("list", n.chunks)
  for (j in 1:n.chunks) {
    cts[[j]] <- cov[tfbs.gr[[j]]]
  }

  # Removing TFBSs that have an outlier coverage
  ct.means <- lapply(cts, mean)
  median.of.ct.means <- median(unlist(ct.means))
  mad.of.ct.means <- mad(unlist(ct.means))

  for (j in 1:n.chunks) {
    outliers <- c(which(ct.means[[j]] < median.of.ct.means - 5 * mad.of.ct.means),
                  which(ct.means[[j]] > median.of.ct.means + 5 * mad.of.ct.means)) 
    if (length(outliers) > 0) {
      cts[[j]] <- cts[[j]][-outliers]
    }
  }

  n.tfbs <- sum(sapply(cts, length))
  rm(ct.means)

  # Computing the relative coverage metric
  mean.cov <- rep(0, 6001)
  for (j in 1:n.chunks) {
    if (length(cts[[j]]) > 0) { # In case all TFBSs from a certain chuk were outliers
      mean.cov <- mean.cov + colSums(as.matrix(cts[[j]]))
    }
  }
  mean.cov <- mean.cov / n.tfbs
  names(mean.cov) <- -3000:3000

  tmp.rel_cov <- ( sum(mean.cov[names(mean.cov) %in% dip_pos]) / length(dip_pos) ) / ( sum(mean.cov[names(mean.cov) %in% ctl_pos]) / length(ctl_pos) )
  rel_cov[i] <- tmp.rel_cov
  rm(mean.cov)
}
names(rel_cov) <- gsub(".rds", "", tfbs.sets)

# Saving the relative coverage metric for each TFBS
saveRDS(rel_cov, file.path(outDir, "tmp", paste0(s, "_rel_cov.rds")))

# Removing big objects
rm(list = c("cov", "cts"))

# If this is the last sample being processed, merge the relative
# coverage for all TFBSs across all samples into a summary table
finished.files <- list.files(file.path(outDir, "tmp"))
n.jobs <- length(frag.files)
n.finished <- length(finished.files)

if (n.jobs == n.finished) {
  # Waiting to make sure the last file is completely written to disk
  Sys.sleep(30)
  
  cov.list <- lapply(finished.files, function(t) readRDS(file.path(file.path(outDir, "tmp", t))))
  cov.df <- do.call("rbind", cov.list)
  rownames(cov.df) <- gsub("_rel_cov.rds$", "", finished.files)

  write.table(cov.df, file = file.path(outDir, "cohort_rel_cov.txt"),
              row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
} else {
  q("no")
}




