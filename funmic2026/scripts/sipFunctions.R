plotPCAWithSpecies <- function (whichPS, scalingType=1, scalingAmount=3, ptitle) {
    sample_data(whichPS)$Isotope = as.factor(sample_data(whichPS)$Isotope)
    ps.hell <- transform_sample_counts(whichPS, function(x) sqrt(x / sum(x)))
    ord.PCA <- ordinate(ps.hell, method="RDA")
    plot.ord <- plot_ordination(ps.hell, ord.PCA, color="Isotope", shape="Isotope")
    speciesMat <- vegan::scores(ord.PCA, display="species", scaling=scalingType)
    labeldf <- data.frame(labels = rownames(speciesMat), speciesMat*scalingAmount)
    label_map <- aes(x = PC1,
        y = PC2,
        shape = NULL,
        color = NULL,
        label = labels,
    )
    plot.ord + geom_point(size=3.5) +
    geom_text( label=sample_data(ps.hell)$Fraction,
        nudge_x = 0.03, nudge_y = 0.03,
        check_overlap = T) + ggtitle(ptitle) +
        geom_text(
          mapping = label_map,
          size = 4,
          data = labeldf,
          show.legend = FALSE,
          check_overlap = T
        )
}

getFractionAbundances = function(ASVname,whichPS,ptitle="",fracOrBD="BD"){
    samp.12 <- rownames(sample_data(whichPS)[sample_data(whichPS)$Isotope == "12",])
    ps.12 <- prune_samples(samp.12, whichPS)
    ps.12.ASV <- prune_taxa(ASVname, ps.12)
    df12 <- as.data.frame(cbind(sample_data(ps.12.ASV)[,fracOrBD], otu_table(ps.12.ASV)))
    df12[] <- sapply(df12, as.numeric)
    colnames(df12) <- c(fracOrBD,"Abundance")
    df12$Isotope <- 12
    ## repeat for 13 isos
    samp.13 <- rownames(sample_data(whichPS)[sample_data(whichPS)$Isotope == "13",])
    ps.13 <- prune_samples(samp.13, whichPS)
    ps.13.ASV <- prune_taxa(ASVname, ps.13)
    df13 <- as.data.frame(cbind(sample_data(ps.13.ASV)[,fracOrBD], otu_table(ps.13.ASV)))
    df13[] <- sapply(df13, as.numeric)
    colnames(df13) <- c(fracOrBD,"Abundance")
    df13$Isotope <- 13
    df12_13 <- rbind(df12,df13)
    df12_13$Isotope <- as.factor(df12_13$Isotope)
    print(tax_table(ps.12.ASV))
    ## how would we do this with ggplot?
    if(fracOrBD == "BD"){
    ggplot(df12_13, aes(x=BD, y=Abundance, color=Isotope, shape=Isotope)) +
      geom_point(size=4) +
      geom_smooth(se=FALSE, fullrange=FALSE, linetype="dashed") + ggtitle(ptitle)
    } else if (fracOrBD == "Fraction"){
    ggplot(df12_13, aes(x=Fraction, y=Abundance, color=Isotope, shape=Isotope)) +
      geom_point(size=4) +
      geom_smooth(se=FALSE, fullrange=FALSE, linetype="dashed") + ggtitle(ptitle)
    }
}


plotFamilies = function(howManyASVs=30){
    topX <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:howManyASVs]
    ps.topX <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
    ps.topX <- prune_taxa(topX, ps.topX)
    plot_bar(ps.topX, x="Fraction", fill="Family") + facet_wrap(~Substrate + Isotope, scales="free_x")
}

NMS_braycurtis = function(physeq){
    sample_data(physeq)$Isotope = as.factor(sample_data(physeq)$Isotope)
    ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
    aa <- plot_ordination(physeq, ord.nmds.bray, color="Substrate", shape="Isotope", title="Bray NMDS")
    aa + geom_point(size=3.5) + geom_text(
                                      label=sample_data(physeq)$Fraction,
                                      nudge_x = 0.05, nudge_y = 0.05,
                                      check_overlap = T)
}
