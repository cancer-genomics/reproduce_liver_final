#!/bin/bash

merg_dir="/users/zfoda/liver_resub/preprocess/tfhi/" 
out_dir="/users/zfoda/liver_resub/preprocess/hiscripts"

cp /dcs04/scharpf/data/zfoda/scripts/v7/submit.sh /users/zfoda/liver_resub/preprocess/hiscripts

mkdir $out_dir/


cd $merg_dir
sample_array=($(ls | sed 's/.*\///g' |sed 's/\.[^.]*$//'))

for sample in ${sample_array[@]}

do
	echo '
library(BSgenome.Hsapiens.UCSC.hg19)


test_args <- function() { index <<- 6}
#test_args()
#Getting command line arguments
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

# Reading in the TFBS region and finding the center of the peak +/- 3kb
liverpeaks.file <- "/users/zfoda/liver_resub/preprocess/tfhi/'${sample}'.rds"
if (!file.exists(liverpeaks.file)) {
  tfbs.gr <- rtracklayer::import("/users/zfoda/liver_resub/preprocess/tfhi/'${sample}'.bed")
  tfbs.gr <- tfbs.gr[seqnames(tfbs.gr) %in% paste0("chr", 1:22)]
  seqlevels(tfbs.gr) <- paste0("chr", 1:22)
  tfbs.gr <- sort(tfbs.gr)
  ranges(tfbs.gr) <- mid(tfbs.gr)
  tfbs.gr <- tfbs.gr + 3000
  names(tfbs.gr) <- tfbs.gr$name
  tfbs.gr <- tfbs.gr[,-1]
  saveRDS(tfbs.gr, file = liverpeaks.file)
} else {
 tfbs.gr <- readRDS(liverpeaks.file)
}




# Reading in the files

fragfiles <- "/dcs04/scharpf/data/annapragada/Liver_Curated/Full_Cohort_Final_Aug22/granges"

fragfile <- list.files(fragfiles, full.names = TRUE)

file <- fragfile[index]
names <- gsub(".rds$", "", basename(file))


# Defining the output file path
ctDir <- "/dcs04/scharpf/data/zfoda/outdir/'${sample}'_cts"
if (!dir.exists(ctDir)) dir.create(ctDir)

filename <- file.path(ctDir, paste0(names, "_cts.rds"))
if(file.exists(filename)) q('"no"')

# Assigning seqlengths to frag.gr
frag.gr <- readRDS(file)
chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
seqinfo(frag.gr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[chrs]


if (!file.exists(filename)) {
  # Computing per base coverage
  cov <- coverage(frag.gr)

  # Creating a matrix to store TTBS counts
  cts <- matrix(nrow = length(tfbs.gr), ncol = length(-3000:3000), dimnames = list(names(tfbs.gr), -3000:3000))
  # Adding per base coverage for each TSS to cts. Reversing the direction of the counts if the gene is on the minus strand.
  for (i in seq_along(tfbs.gr)) {
    print(i)
    try({bp.cov <- as.numeric(cov[[as.character(seqnames(tfbs.gr[i]))]][start(tfbs.gr[i]):end(tfbs.gr[i])])
    cts[i,] <- bp.cov
    print(paste0(i," good"),quote = FALSE)
    })
  }
  saveRDS(cts, file = filename)
}
' >$out_dir/${sample}_coverage-03.R

echo '
#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=16G
#$ -l h_vmem=16G
#$ -t 1-560
#$ -N '${sample}'
#$ -o ../logs2/


module load conda_R/4.0.x

Rscript ./'${sample}'_coverage-03.R $SGE_TASK_ID
' >$out_dir/${sample}_coverage-03.sh


echo '
# Request 96GB of RAM

library(dplyr)
library(ggplot2)
library(pROC)
library(readr)
library(readxl)

if (!file.exists("/users/zfoda/liver_resub/preprocess/outdirhigh/roc'${sample}'.pdf")){
# Reading in the metadata
meta<-read_excel("/users/zfoda/liver_resub/data/Clinical_Metadata_spreadsheet_8_11.xlsx") #%>% mutate(`Sex(1=F,2=M)` = replace(`Sex(1=F,2=M)`, `Sex(1=F,2=M)` =="1", "Female"))%>% mutate(`Sex(1=F,2=M)` = replace(`Sex(1=F,2=M)`, `Sex(1=F,2=M)` =="2", "Male")) %>% rename("Sex" =`Sex(1=F,2=M)`,"id"=CGID)

 
#finding cts files
ctDir <- "/dcs04/scharpf/data/zfoda/outdir/'${sample}'_cts"
ct.files <- list.files(ctDir)
ids <-  gsub("_cts.rds", "", ct.files)
#ids <- gsub("CGPLH921P_1", "CGPLH921P", ids)
#loading and combining cts

df <- data.frame(id = character(0), grp = character(0), pos = character(0), cov= character(0), rel_cov = character(0))



for (i in 1:length(ct.files)) {
  try({rds.tmp <- readRDS(file.path(ctDir, ct.files[i]))
  ct <- colMeans(rds.tmp, na.rm = TRUE)
  tmp.df <- data.frame(id = ids[i],
                       grp = meta$`Disease`[meta$id == ids[i]],
                       pos = as.numeric(names(ct)),
                       cov = as.numeric(ct),
                       rel_cov = as.numeric(ct) / max(as.numeric(ct)))
  df <- rbind(df, tmp.df)
  rm(rds.tmp)
  })
}

# Recoding grp as HCC or Other
df$grp1 <- ifelse(df$grp == "HCC", "HCC", "Other")

# Saving the raw data
write.table(df, file = "/dcs04/scharpf/data/zfoda/outdir_resub/'${sample}'_coverage1.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


# Computing relative coverage
ids <- unique(df$id)
df$rel_cov <- NA
for (i in 1:length(ids)) {
  id.hits <- which(df$id == ids[i])
  df$rel_cov[id.hits] <- df$cov[id.hits] / max(df$cov[id.hits])
}

# Plotting relative coverage for all samples separately
df$grp1 <- factor(df$grp1, levels = c("Other", "HCC"))
id_ord <- c(unique(subset(df, grp1 == "Other")$id), unique(subset(df, grp1 == "HCC")$id))
df$id <- factor(df$id, levels = id_ord)

# Plotting the non-HCC group as the median and the 0.05 and 0.95 quantiles
q.df <- df %>% 
          filter(grp1 != "HCC") %>% 
          group_by(pos) %>% 
          summarise(q5 = quantile(rel_cov, 0.05),
                    m = median(rel_cov),
                    q95 = quantile(rel_cov, 0.95)) %>%
          as.data.frame

hcc.df <- subset(df, grp1 == "HCC")

# Adding the median of the 'Other' group to sclc.df
median.df <- data.frame(id = "Other_median",
                        grp= "Other",
                        grp1 = "Other",
                        pos = q.df$pos,
                        cov = NA,
                        rel_cov = q.df$m)
hcc.df <- rbind(median.df, hcc.df)
hcc.df$grp1 <- factor(hcc.df$grp1, levels = c("HCC", "Other"))

# Plotting the relative coverage at the sites
#p <- 	hcc.df  %>% ggplot(aes(x = pos, y = rel_cov, color = grp, group = id)) +
#       geom_ribbon(data = q.df, aes(x = pos, ymin = q5, ymax = q95), inherit.aes = FALSE, alpha = 0.5, color = NA, fill = "grey50") +
 #      geom_line(aes(color = grp), alpha = 0.8) +
  #     cowplot::theme_cowplot() +
   #    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  #     ylim(c(0.5, 1)) +
  #     ylab("Relative coverage") +
  #     xlab("Position relative to peak (bp)") +
  #     scale_color_manual(values = c(HCC = "slateblue", Other = "black")) +
  #     theme(legend.title = element_blank(),
  #           legend.position = c(0.8, 0.2))

#ggsave(p, filename = "/users/zfoda/liver_resub/preprocess/outdirhigh/'${sample}'_coverage1.pdf",
   #    height = 4, width = 8)

# Computing a relative coverage metric for each sample for the ROC curve
ctl_pos <- c(-3000:-2500, 2500:3000)
dip_pos <- -200:200
summary.df <- data.frame(id = character(0),
                         grp = character(0),
                         grp1 = character(0),
                         rel_cov = character(0),
                         TF =character(0))
for (i in 1:length(ids)) {
  z <- subset(df, id == ids[i])
  flank_height <-((sum(head(sort(z$cov[z$pos %in% dip_pos], decreasing = TRUE),50)))/50)
  background <- (sum(z$cov[z$pos %in% ctl_pos]) / length(ctl_pos))
  base_cov <- (sum(tail(sort(z$cov[z$pos %in% dip_pos], decreasing = TRUE),10)))/10
  access <- log2(flank_height/background)
  depth <- log2(base_cov/flank_height)
  rel_cov <- ( sum(z$cov[z$pos %in% dip_pos]) / length(dip_pos) ) / ( sum(z$cov[z$pos %in% ctl_pos]) / length(ctl_pos) )  
  tmp.summary.df <- data.frame(id = ids[i],
                               grp = unique(df$grp[df$id == ids[i]]),
                               grp1 = unique(df$grp1[df$id == ids[i]]),
                               rel_cov = rel_cov,
                               access =access,
                               depth= depth,
                               TF = "'${sample}'")
  summary.df <- rbind(summary.df, tmp.summary.df)
}

write.table(summary.df, file = "/users/zfoda/liver_resub/preprocess/outdirhigh/'${sample}'_summary.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Creating the ROC curve
roc.res <- roc(summary.df$grp1, summary.df$rel_cov, ci = TRUE)
roc.res1 <- roc(summary.df$grp1, summary.df$depth, ci = TRUE)
pdf(file = "/users/zfoda/liver_resub/preprocess/outdirhigh/roc'${sample}'.pdf", width = 5, height = 5)
plot(roc.res, print.auc = TRUE, print.auc.x = 0.5, print.auc.y = 0.35, print.auc.pattern = "CovAUC: %.2f (%.2f-%.2f)")
plot(roc.res1, print.auc = TRUE, print.auc.x = 0.5, print.auc.y = 0.25, print.auc.pattern = "DepthAUC: %.2f (%.2f-%.2f)")

dev.off()

}' >$out_dir/plot_${sample}_coverage-04.R

echo '
#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_fsize=1000G
#$ -o ../logs2/


module load conda_R/4.0.x

Rscript ./plot_'${sample}'_coverage-04.R' >$out_dir/plot_${sample}_coverage-04.sh


done

