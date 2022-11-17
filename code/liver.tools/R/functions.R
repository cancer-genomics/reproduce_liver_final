#' @include help.R
NULL

###########
#FXNs
.meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
    meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
        n <- length(x)
        y <- rep(NA,n)

        window.mean <- function(x, j, k, na.rm=na.rm){
            if (k>=1){
                return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
            } else {
                return(x[j])
            }
        }

        for (i in (k+1):(n-k)){
            y[i] <- window.mean(x,i,k, na.rm)
        }
        for (i in 1:k){
            y[i] <- window.mean(x,i,i-1, na.rm)
        }
        for (i in (n-k+1):n){
            y[i] <- window.mean(x,i,n-i,na.rm)
        }
        y
    }

    for (i in 1:iter){
        x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
    }
    x
}

geom_ab_sc_lymph_agree <- function(tib, chromosome) {
    tib  <- tib %>% filter(chr==chromosome) %>%
         mutate(eigen = .meanSmoother(eigen),
                color=ifelse(eigen < 0, "gray50", "brickred"))
    p <- ggplot(tib, aes(x=bin, y=eigen,
                         color=lymph_agree,fill=lymph_agree)) +
        geom_bar(stat="Identity", width=1)
    p <- p + facet_wrap(vars(source2), ncol=1)
    p <- p +  scale_x_continuous(expand = c(0, 0)) +
        xlab(chromosome) +theme_classic()+
        theme(legend.position="none",axis.title.y=element_blank()) +
        scale_color_manual(values=c("red4", "gray50"))+
        scale_fill_manual(values=c("red4", "gray50"))
    p
}

geom_ab <- function(tib, chromosome) {
    tib  <- tib %>% filter(chr==chromosome) %>%
         mutate(eigen = .meanSmoother(eigen),
                color=ifelse(eigen < 0, "gray50", "brickred"))
    p <- ggplot(tib, aes(x=bin, y=eigen, color=color,fill=color)) +
        geom_bar(stat="Identity", width=1)
    p <- p + facet_wrap(vars(source2), ncol=1)
    p <- p +  scale_x_continuous(expand = c(0, 0)) +
        xlab(chromosome) +theme_classic() +
        theme(legend.position="none",axis.title.y=element_blank()) +
        scale_color_manual(values=c("red4", "gray50"))+
        scale_fill_manual(values=c("red4", "gray50"))
    p
}

geom_ab_transp <- function(tib, chromosome,mis_bins) {
    tib  <- tib %>% filter(chr==chromosome) %>%
        mutate(color=ifelse(eigen < 0, "gray50", "red4"),
               transp=ifelse(bin %in% mis_bins,.95,.1))
    p <- ggplot(tib, aes(x=bin, y=eigen,fill=color,alpha=transp)) +
        geom_bar(stat="Identity", width=1)
    p <- p + facet_wrap(vars(source2), ncol=1)
    p <- p +  scale_x_continuous(expand = c(0, 0)) +
        xlab(chromosome) +theme_classic()+
        theme(legend.position="none",axis.title.y=element_blank()) +
        scale_color_manual(values=c("red4", "gray50"))+
        scale_fill_manual(values=c("red4", "gray50")) +
        scale_alpha_identity()
    p
}

geom_ab_transp_overlay <- function(tib, chromosome,mis_bins) {
    tib  <- tib %>% filter(chr==chromosome) %>%
        mutate(color=ifelse(eigen < 0, "gray50",
                            "red4"),transp=ifelse(bin %in%
                                                  mis_bins,.95,.1))
    p <- ggplot(tib,
                aes(x=bin, y=eigen,color=source2,
                    fill=color,alpha=transp)) +
        geom_bar(stat="Identity", width=1)
    ##p <- p + facet_wrap(vars(source2), ncol=1)
    p <- p +  scale_x_continuous(expand = c(0, 0)) +
        xlab(chromosome) +theme_classic()+
        theme(legend.position="none",axis.title.y=element_blank()) +
        scale_color_manual(values=c("green","black","red4", "gray50"))+
        scale_fill_manual(values=c("red4", "gray50")) +
        scale_alpha_identity()
    p
}

#' @export
chr22_wrangling <- function(dat){
    sources <- unique(dat$source2) %>% sort() %>%
        "["(!grepl("[357]", .))
    slevels <- sources[c(1, 2, 4, 3, 5)]
    ## TODO: replace "Healthy" with "Noncancer"
    labels <- c("HCC Reference AB Compartments",
                "Extracted HCC Plasma Component",
                "HCC Plasma",
                "Healthy Plasma",
                "Lymphoblast HiC Reference AB Compartments")
    c16 <- dat %>% filter(chr=="chr22") %>%
        filter(source2 %in% sources) %>%
        mutate(source2=factor(source2, slevels, labels))
    ##the bins that don't agree b/w reference tracks
    mis_bins <- dat %>%
        filter(source2=="1. AB Compartments LIHC" & lymph_agree=="no")
    mis_bins <- mis_bins$bin
    ## trying to do victor's thing of bins with a big change in magnitude
    lymph <- dat %>% filter(source2=="8. EBV transformed lymhoblast HiC")
    lihc <- dat %>% filter(source2=="1. AB Compartments LIHC")
    bins <- inner_join(lymph, lihc, by="bin")
    bins$diff <- bins$eigen.x - bins$eigen.y
    bins$z_score <- (bins$diff-mean(bins$diff))/sd(bins$diff)
    sig_bins <- bins %>% filter(z_score > 1.96 | z_score< -1.96)
    ##p=.05, maybe should be corrected but not really ind tests?
    sig_bins <- sig_bins$bin
    d <- setdiff(sig_bins, mis_bins)
    ## adding the discordant sign bins to the track and make the informative bins stand out
    d <- append(d, mis_bins)
    ##trace(geom_ab_transp, browser)
    ##b <- geom_ab_transp(c16, "chr22", d)
    tib <- c16
    chromosome <- "chr22"
    mis_bins <- d
    track.data <- tib %>% filter(chr==chromosome) %>%
        mutate(color=ifelse(eigen < 0, "gray50", "red4"),
               transp=ifelse(bin %in% mis_bins, .95, .1))
    ##now making the other panel thing to show mag diff
    c16 <- dat %>% filter(chr=="chr22")
    ##only the plasma tracks
    c16 <- c16 %>%
        filter(source2 %in% c("4. Median Healthy Plasma Ratios",
                              "6. Median Ichor High Plasma Ratios"))
    c16$source2 <- factor(c16$source2,
                          levels=c("6. Median Ichor High Plasma Ratios",
                                   "4. Median Healthy Plasma Ratios"),
                          labels=c("HCC Plasma","Healthy Plasma"))
                                        #coloring the bins as above
    tib <- c16 %>% filter(chr=="chr22") %>%
        mutate(eigen = .meanSmoother(eigen),
               color=ifelse(eigen < 0, "gray50", "red4"),
               transp=ifelse(bin %in% d,.95,.1))
    ##the title thing is weird and it needs a legend also the lines
    ##are too thick, fix later
    tib$source3 <- "Plasma -- White = HCC, Black = Healthy"
    list(tib=tib, track.data=track.data)
}

#' @export
deconvolution_boxplots <- function(data){
    mis_bins <- data %>%
        filter(source2=="1. AB Compartments LIHC" & lymph_agree=="no")
    mis_bins <- mis_bins$bin
    data2 <- data %>% filter(bin %in% mis_bins)
    ##adding a bunch of columns that could be cool
    lymph_agreement <- data2 %>% group_by(source2,chr) %>%
        count(lymph_agree=="yes") %>%
        filter(`lymph_agree == "yes"`==TRUE)
    lihc_agreement <- data2 %>%
        group_by(source2,chr) %>%
        count(lihc_agree=="yes") %>%
        filter(`lihc_agree == "yes"`==TRUE)
    all_groups <- data2 %>%
        group_by(source2,chr) %>% count()
    lymph_agreement$combo <- paste(lymph_agreement$source2,
                                   lymph_agreement$chr,sep=" ")
    lihc_agreement$combo <- paste(lihc_agreement$source2,
                                  lihc_agreement$chr,sep=" ")
    all_groups$combo <- paste(all_groups$source2,
                              all_groups$chr, sep=" ")
    lymph_agreement <- inner_join(lymph_agreement,
                                  all_groups,by="combo")
    lymph_agreement$percentage <- lymph_agreement$n.x/lymph_agreement$n.y
    lymph_agreement$type <- "Lymphoblast Agreement"
    lymph_agreement <- lymph_agreement %>%
        select(source2.x,combo,chr.x,percentage,type)
    lihc_agreement <- inner_join(lihc_agreement, all_groups,by="combo")
    lihc_agreement$percentage <- lihc_agreement$n.x/lihc_agreement$n.y
    lihc_agreement$type <- "LIHC Agreement"
    lihc_agreement <- lihc_agreement %>%
        select(source2.x,combo,chr.x,percentage,type)
    agree <- rbind(lihc_agreement,lymph_agreement)
    lihc_agreement <- lihc_agreement %>%
        rename(Liver_Percentage=percentage)
    lymph_agreement <- lymph_agreement %>%
        rename(Lymphoblast_Percentage=percentage)
    perc <- inner_join(lihc_agreement %>%
                       select(Liver_Percentage, combo, source2.x, chr.x),
                       lymph_agreement %>%
                       select(Lymphoblast_Percentage,combo),by="combo")
    perc2 <- perc %>%
        filter(source2.x=="2. Extracted Liver component - Ratios" |
               source2.x=="4. Median Healthy Plasma Ratios" |
               source2.x=="6. Median Ichor High Plasma Ratios" )
    chrlevels <- paste0("chr", 1:22)
    perc2$chr.x <- factor(perc2$chr.x, chrlevels)
    slevels <- c("2. Extracted Liver component - Ratios",
                 "6. Median Ichor High Plasma Ratios",
                 "4. Median Healthy Plasma Ratios")
    perc2 <- perc2 %>%
        mutate(source2.x=factor(source2.x,
                                slevels,
                                c("Estimated Liver Component",
                                  "HCC Plasma",
                                  "Healthy Plasma")),
               odds=(Liver_Percentage)/(1-Liver_Percentage))
    perc2
}

#' @export
zscore_features <- function(features, meta){
    features <- features %>% select(id,starts_with("z"))
    features <- inner_join(features, meta %>%
                                     select(id,`HCC Status`,Disease),
                           by=c("id"="id"))
    features <- features %>% filter(Disease != "Cholangiocarcinoma")
    features <- features %>%
        gather(key="Arm", value="Zscore", starts_with("z"))
    features$Arm <- sapply( str_split(features$Arm,"_"), "[", 2 )

    chrlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p",
                   "5q","6p","6q","7p","7q","8p","8q","9p","9q",
                   "10p","10q","11p","11q","12p","12q",
                   "13q","14q","15q","16p","16q","17p","17q",
                   "18p","18q","19p","19q","20p","20q","21q","22q")
    features$Arm <- factor(features$Arm,
                           levels = chrlevels)
    features <- features %>% mutate(colors=if_else(Zscore<=0,"Negative","Positive"))
    ##features$colors<-factor(features$colors,levels=c("Positive","Negative"))
    ##Z-score for two-tailed p-value .05 w/ bonferroni correction for 39 arms
    features <- features %>%
        mutate(colors=if_else(Zscore>= -3.22 & Zscore <= 3.22 ,
                              "Neutral", features$colors))
    ##features<- features %>% mutate(transp=if_else(colors=="Neutral",.9,.2))
    ##features<- features %>% mutate(transp=if_else(colors=="Neutral",.05,.99))
    features <- features %>% mutate(transp=if_else(colors=="Neutral",1,1))
    ##features<-features %>% mutate(colors=if_else(Zscore<=0,"Negative","Positive"))
    features$colors<-factor(features$colors,levels=c("Positive","Negative","Neutral"))
    features$y=1
    ##features$log<-log10(abs(features$Zscore))
    ##features <- features %>% mutate(log=if_else(Zscore<=0,-1*log,log))
    features <- features %>% mutate(Zscore=if_else(Zscore<=-150,-150,Zscore))
    features <- features %>% mutate(Zscore=if_else(Zscore>150,150,Zscore))
    features <- features %>%
        mutate(`HCC Status`=if_else(`HCC Status`=="Yes","HCC","Viral Hepatitis\n and Healthy"))
    features <- features %>%
        mutate(`HCC Status`=if_else(Disease=="Cirrhosis","Cirrhosis",`HCC Status`))
    features$hcc_status <- factor(features$`HCC Status`,
                                  levels=c("HCC","Cirrhosis","Viral Hepatitis\n and Healthy"))
    features2  <- features %>%
        rename(z=Zscore) %>%
        mutate(root2=nroot(z, 2)) %>%
        mutate(root3=nroot(z, 3))
    features2
}

#' @export
top_annotation <- function(sc){
    ta <- subset(sc, grepl('PC', features)) %>% filter(feature.type=="ratio")
    ta$pc.id <- gsub('PC', '', gsub(' ','',ta$features))
    ta$sign.value = factor(ta$sign.value, levels = c(-1, 1))
    ta <- ta[,c('pc.id', 'abs.value', 'sign.value')]
    #ta <- rbind(ta, data.frame(pc.id = setdiff(seq(3), ta$pc.id), abs.value =0, sign.value = NA))
    ta$pc.id <- factor(ta$pc.id, levels = ta$pc.id)
    ta$x = 1
    ta$arm = 1
    ta
}


#' @export
side_annotation <- function(pd, sc, bins5mb){
    ddply <- plyr::ddply
    . <- plyr::.
    locs <- unique(bins5mb[,c('bin', 'chr', 'start', 'end', 'arm')])
    locs$pos <- apply(locs[,c('start', 'end')], 1, mean)
    bl <- unique(pd[,c('chr','pos', 'bin.id', 'arm')])
    bl$x = 1
    arm.imp <- ddply(bl, .(arm), summarize,
                     bin.id = mean(bin.id), x = mean(x))
    arm.coef <- subset(sc, ! grepl('PC', features)) %>% filter(feature.type=="zscore") %>%
        "["(,c('features', 'abs.value', 'sign.value')) 
    arm.coef$arm <- gsub('Z ', '', arm.coef$features)
    arm.coef$features <- NULL
    arm.coef <- rbind(arm.coef,
                      data.frame(abs.value = 0, sign.value = NA,
                                 arm = setdiff(unique(locs$arm),
                                               arm.coef$arm)))
    arm.imp <- merge(arm.imp,
                     arm.coef[,c('arm', 'abs.value', 'sign.value')],
                     by = 'arm', all.x = TRUE)
    bl$pc.id <- 1
    list(bl=bl,
         arm.imp=arm.imp)
}

#' @export
side_annotation <- function(pd, sc, bins5mb){
    ddply <- plyr::ddply
    . <- plyr::.
    locs <- unique(bins5mb[,c('bin', 'chr', 'start', 'end', 'arm')])
    locs$pos <- apply(locs[,c('start', 'end')], 1, mean)
    bl <- unique(pd[,c('chr','pos', 'bin.id', 'arm')])
    bl$x = 1
    arm.imp <- ddply(bl, .(arm), summarize,
                     bin.id = mean(bin.id), x = mean(x))
    arm.coef <- subset(sc, ! grepl('PC', features)) %>% filter(feature.type=="zscore") %>%
        "["(,c('features', 'abs.value', 'sign.value')) 
    arm.coef$arm <- gsub('Z ', '', arm.coef$features)
    arm.coef$features <- NULL
    arm.coef <- rbind(arm.coef,
                      data.frame(abs.value = 0, sign.value = NA,
                                 arm = setdiff(unique(locs$arm),
                                               arm.coef$arm)))
    arm.imp <- merge(arm.imp,
                     arm.coef[,c('arm', 'abs.value', 'sign.value')],
                     by = 'arm', all.x = TRUE)
    bl$pc.id <- 1
    list(bl=bl,
         arm.imp=arm.imp)
}

#' @export
tf_annotation <- function(pd, sc, bins5mb){
    ddply <- plyr::ddply
    . <- plyr::.
    locs <- unique(bins5mb[,c('bin', 'chr', 'start', 'end', 'arm')])
    locs$pos <- apply(locs[,c('start', 'end')], 1, mean)
    bl <- unique(pd[,c('chr','pos', 'bin.id', 'arm')])
    bl$x = 1
    arm.imp <- ddply(bl, .(arm), summarize,
                     bin.id = mean(bin.id), x = mean(x))
    arm.coef <- sc %>% filter(feature.type=="TFBS") %>%
        "["(,c('features', 'abs.value', 'sign.value')) 
    arm.coef$arm <- sapply(str_split(arm.coef$features,"[.]"),"[",2)
    arm.coef$features <- NULL
    arm.coef$x<-1
    arm.coef$bin.id<-c(130,140,150,160,170,180,190,200,210)
    bl$pc.id <- 1
    list(bl=bl,
         arm.imp=arm.coef)
}


nroot <- function(x, n){
    e <- 1/n
    posx <- x > 0
    negx <- x < 0
    x[posx] <- x[posx]^e
    x[negx] <- -(-x[negx])^e
    x
}

random_sample <- function(iter, x, n){
    x <- filter(x, hcc_status == "HBV and Healthy")
    id_ <- sample(unique(x$id), n, replace=FALSE)
    x %>%
        filter(id %in% id_) %>%
        arrange(Arm, root2) %>%
        mutate(iter=iter) %>%
        select(hcc_status, Arm, root2, root3, iter)
}

highlight <- function(x){
    xx <- x %>%
        mutate(hcc_z=root3^3,
               nc_z=h_r3^3,
               delta=hcc_z - nc_z) %>%
        mutate(signif=ifelse(delta < -5, "Deletion", "Neutral"),
               signif=ifelse(delta > 5, "Amplification", signif)) %>%
        ##mutate(d=root2 - h_r2,
        ##       is_2sd=abs(d) > 2,
        ##       direction=sign(root2),
        ##       signif=ifelse(is_2sd & direction==1,
        ##                     "Amplification", "Neutral"),
        ##       signif=ifelse(is_2sd & direction==-1,
        ##                     "Deletion", signif)) %>%
        pull(signif)
    xx
}

## params: named list with following elements
##  - prevalence
##  - sens (sensitivity)
##  - spec (specificity)
##  - compliance
performance <- function(prevalence, adherence, tool, Sensitivity,
                        Specificity){
    sens <- Sensitivity
    spec <- Specificity
    L <- length(adherence)
    screened <- rbinom(L, size=100e3, prob=adherence)
    ## prevalence
    P <- rbinom(L, size=screened, prob=prevalence) ## P = FN + TP
    N <- screened - P  ## N = TN + FP
    TP <- rbinom(L, size=P, prob=sens)
    FP <- rbinom(L, size=N, prob=(1-spec))
    TN <- N - FP  ## N = FP + TN
    FN <- P - TP
    fpr <- FP/N
    fnr <- FN/P
    tnr <- TN/N
    tpr <- TP/P
    acc <- (TP+TN)/(P+N)
    err <- (FP+FN)/(P+N)
    ppv <- TP/(TP+FP)
    npv <- TN/(TN+FN)
    stats <- tibble("P"=P, "N"=N, "TP"=TP, "FP"=FP, "TN"=TN,
          "FN"=FN, "acc"=acc, "err"=err,
          "fpr"=fpr, "fnr"=fnr, "tnr"=tnr,
          "tpr"=tpr, "ppv"=ppv, "npv"=npv,
          "number_screened"=screened,
          device=tool)
    stats
}
