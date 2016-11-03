# Bar plot of community structure

#Colin Brislawn, 2015-12-22
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer", "scales")
library("vegan")

theme_set(theme_bw())

#####
# Importing data

no_meta <- import_biom(file.path('../../data/otus_vsearch/otu_table_w_tax.biom'),
											 file.path('../../data/otus_vsearch/rep_set.tre'))
no_meta

meta <- import_qiime_sample_data(file.path('../../metadata/metadata.txt'))

full <- merge_phyloseq(meta, no_meta)
full



#####
# Identify the important variables
metadata <- sample_data(full)
head(metadata)
metadata.n_unique <- rapply(metadata, function(x) length(table(x)))
# The pairwise interesting factors
metadata.n_unique[metadata.n_unique == 2]
# Oh yeah, there are no categorical variables in this study. All continuous, all the time. 
# Other potentially interesting factors 
metadata.n_unique[metadata.n_unique > 3 & metadata.n_unique < max(metadata.n_unique)]

# So interesting factors are 
#MatDepth.cryosection.rel    Gross.PS.rate 
#pctC.replaced  pctC.replaced.perhr  X13C.isotoperatio.rate  irradiance.PAR    dissolved.ox  porosity 


######
# Tell R to treat numbers as numbers (not factors)
sample_data(full)$MatDepth.cryosection.mm <- as.numeric(levels(sample_data(full)$MatDepth.cryosection.mm))[sample_data(full)$MatDepth.cryosection.mm]
sample_data(full)$Gross.PS.rate <- as.numeric(levels(sample_data(full)$Gross.PS.rate))[sample_data(full)$Gross.PS.rate]
sample_data(full)$MatDepth.cryosection.rel <- as.numeric(levels(sample_data(full)$MatDepth.cryosection.rel))[sample_data(full)$MatDepth.cryosection.rel]

sample_data(full)$pctC.replaced <- as.numeric(levels(sample_data(full)$pctC.replaced))[sample_data(full)$pctC.replaced]
sample_data(full)$pctC.replaced.perhr <- as.numeric(levels(sample_data(full)$pctC.replaced.perhr))[sample_data(full)$pctC.replaced.perhr]
sample_data(full)$X13C.isotoperatio.rate <- as.numeric(levels(sample_data(full)$X13C.isotoperatio.rate))[sample_data(full)$X13C.isotoperatio.rate]
sample_data(full)$irradiance.PAR <- as.numeric(levels(sample_data(full)$irradiance.PAR))[sample_data(full)$irradiance.PAR]
sample_data(full)$dissolved.ox <- as.numeric(levels(sample_data(full)$dissolved.ox))[sample_data(full)$dissolved.ox]
sample_data(full)$porosity <- as.numeric(levels(sample_data(full)$porosity))[sample_data(full)$porosity]

head(sample_data(full))


#####
# Filtering and normalization

head(sort(sample_sums(full)))
tail(sort(sample_sums(full)))
# Good depth!

# remove failed samples.
filtered <- prune_samples(sample_sums(full)>=22000, full)
nsamples(full) - nsamples(filtered)
head(sample_data(filtered))
# Like we said, all samples have good depth. We only lose one.

# We will not rarify (we will scale by relative abundance after glomming instead)

#####
# Go to folder for saving graphs
setwd("./barplot/")
list.files()





###########
# Barplot using binned depths.

#####
# Select one cohort
#cohort.rarefied <- subset_samples(filtered.rarefied)                      # All
#cohort.rarefied <- subset_samples(filtered, COHORT %in% c("12A", "12B"))  
cohort <- subset_samples(filtered, !is.na(MatDepth.cryosection.rel))       # Mat samples
cohort

# Make a new column with rounded depth
sample_data(cohort)$MatDepth.cryosection.rel.round <- round(sample_data(cohort)$MatDepth.cryosection.rel, digits = 1)
# Convert depth columsn to factors (for graphing)
#sample_data(cohort)$MatDepth.cryosection.rel.round <- as.factor(sample_data(cohort)$MatDepth.cryosection.rel.round)


#top <- prune_taxa(names(sort(taxa_sums(cohort.rarefied), TRUE))[1:50], cohort.rarefied)
top <- cohort #All taxa! (Good if you want to glom before getting the top ones)
top

#merge samples by a category to increase legibility
#top.merge <- merge_samples(top, group = "LinkerPrimerSequence") # This merges them into one.
top.merge <- merge_samples(top, group = "MatDepth.cryosection.rel.round") # Merge samples that appear at the same depth
#top.merge <- top #do not merge
top.merge
head(sample_data(top.merge))


#glom taxa can simplify the graph (but changes the meaning of each internal bar)
top.merge.glom <- tax_glom(top.merge, taxrank = "Rank3")
#top.merge.glom <- top.merge
top.merge.glom
head(tax_table(top.merge.glom))

# Transform to Relative Abundance so that the bars are the same height. 
top.merge.glom <- transform_sample_counts(top.merge.glom, function(x) x / sum(x) )


#we can filter again, after glomming.
#top.merge.glom.top <- top.merge.glom
top.merge.glom.top <- prune_taxa(names(sort(taxa_sums(top.merge.glom), TRUE))[0:11], top.merge.glom)
sum(taxa_sums(top.merge.glom.top)) / sum(taxa_sums(top.merge.glom))

#top.merge.glom.top <- top.merge.glom
top.merge.glom.top
#tax_table(top.merge.glom.top)
#otu_table(top.merge.glom.top)


plot <- plot_bar(top.merge.glom.top, x = "MatDepth.cryosection.rel.round", fill = "Rank3") +
	scale_fill_manual(values = brewer.pal(11, "Spectral"))
plot


# Better graphing
physeqdf <- psmelt(top.merge.glom.top)
physeqdf$glom_rank <- factor(physeqdf[["Rank3"]])
p <- ggplot(data=physeqdf[order(physeqdf$glom_rank),], aes(x=MatDepth.cryosection.rel.round, y=Abundance, fill=Rank3))
p <- p + geom_bar(stat="identity", colour="black", position = "stack")
p


p.bar <- p + scale_fill_brewer(type = 'seq', palette = "Spectral", name = "Taxa") + 
	ylab("Relative Abundance of Common \n Taxa at the Class Level") +
#	theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
	coord_flip() + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1.0), trans="reverse") +
	xlab("Relative Depth")
p.bar

# one column figure
ggsave('Classes by rounded relative depth rotate.pdf', p.bar, width = 86, height = 60, units = 'mm', dpi = 300, scale = 2)
# two columsn figures are 178 mm
# max depth is 210 mm








#"Let’s kick out a bar plot that shows abundance at the phylum level, same format and assumptions made for the current one made at the class level."
# Yes, let's.

top.merge.glom <- tax_glom(top.merge, taxrank = "Rank2")
head(tax_table(top.merge.glom))

# Transform to Relative Abundance so that the bars are the same height. 
top.merge.glom <- transform_sample_counts(top.merge.glom, function(x) x / sum(x) )

#we can filter after glomming.
top.merge.glom.top <- prune_taxa(names(sort(taxa_sums(top.merge.glom), TRUE))[0:11], top.merge.glom)
sum(taxa_sums(top.merge.glom.top)) / sum(taxa_sums(top.merge.glom))

top.merge.glom.top
#tax_table(top.merge.glom.top)
#otu_table(top.merge.glom.top)

# Better graphing
physeqdf <- psmelt(top.merge.glom.top)
physeqdf$glom_rank <- factor(physeqdf[["Rank2"]])
physeqdf$glom_rank
# Trying to reorder stacked bar plot such that taxa are in the same order in every column
# instead of sorting them by abundance and using a different order like geom_bar wants to do.
factor(physeqdf$glom_rank, levels = names(sort(physeqdf$glom_rank)))
head(physeqdf)

p <- ggplot(physeqdf)
p <- ggplot(data=physeqdf[order(physeqdf$glom_rank),], aes(x=MatDepth.cryosection.rel.round, y=Abundance, fill=Rank2))
p <- p + geom_bar(stat="identity", colour="black", position = "stack")
p
#p + facet_grid(CD_RESECTION + PATIENT_NUMBER + IBD_TYPE~.) + scale_fill_manual(values = repeatPallette(n_unique_taxa))

p.bar <- p + scale_fill_brewer(type = 'seq', palette = "Spectral", name = "Taxa") + 
	ylab("Relative Abundance of Common \n Taxa at the Phylum level") +
	#	theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
	coord_flip() + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1.0), trans="reverse") +
	xlab("Relative Depth")
p.bar

# one column figure
ggsave('Phylum by rounded relative depth rotate.pdf', p.bar, width = 86, height = 60, units = 'mm', dpi = 300, scale = 2)
# two columns figures are 178 mm
# max depth is 210 mm







#####
# Higher resolution bar plot so we can investigate more specific taxa
top.merge.glom <- tax_glom(top.merge, taxrank = "Rank6")
head(tax_table(top.merge.glom))

# Transform to Relative Abundance so that the bars are the same height. 
top.merge.glom <- transform_sample_counts(top.merge.glom, function(x) x / sum(x) )

#we can filter after glomming.
top.merge.glom.top <- prune_taxa(names(sort(taxa_sums(top.merge.glom), TRUE))[0:50], top.merge.glom)
sum(taxa_sums(top.merge.glom.top)) / sum(taxa_sums(top.merge.glom))

top.merge.glom.top
#tax_table(top.merge.glom.top)
#otu_table(top.merge.glom.top)

# Better graphing
physeqdf <- psmelt(top.merge.glom.top)
physeqdf$glom_rank <- factor(physeqdf[["Rank2"]])
physeqdf$glom_rank
# Trying to reorder stacked bar plot such that taxa are in the same order in every column
# instead of sorting them by abundance and using a different order like geom_bar wants to do.
factor(physeqdf$glom_rank, levels = names(sort(physeqdf$glom_rank)))
head(physeqdf)

p <- ggplot(data=physeqdf[order(physeqdf$glom_rank),], aes(x=MatDepth.cryosection.rel.round, y=Abundance, fill=Rank2))
p <- p + geom_bar(stat="identity", colour="black", position = "stack")
p
#p + facet_grid(CD_RESECTION + PATIENT_NUMBER + IBD_TYPE~.) + scale_fill_manual(values = repeatPallette(n_unique_taxa))

p.bar <- p + scale_fill_brewer(type = 'seq', palette = "Spectral", name = "Taxa") + 
	ylab("Relative Abundance of Common \n Taxa at the Phylum level") +
	#	theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
	coord_flip() + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1.0), trans="reverse") +
	xlab("Relative Depth")
p.bar

ggsave('highres genus color phylum.pdf', p.bar, width = 200, height = 60, units = 'mm', dpi = 300, scale = 2)
# Also export table
output <- physeqdf[order(physeqdf$Sample, physeqdf$glom_rank),]
head(output)
write.table(output, "hires genus.txt", F, F, "\t", row.names = F)
#










#######
# Export flat OTU tables in tsv format at various depths
top.merge
# The table in question is the one used in the barplot.
# Samples have been merged to the nearest depth of 10% for a 
# total of 11 samples (0, .1, .2... .9, 1.0)

top.merge.Rank2 <- tax_glom(top.merge, taxrank = "Rank2")
top.merge.Rank3 <- tax_glom(top.merge, taxrank = "Rank3")
top.merge.Rank4 <- tax_glom(top.merge, taxrank = "Rank4")
top.merge.Rank5 <- tax_glom(top.merge, taxrank = "Rank5")
#top.merge.Rank6 <- tax_glom(top.merge, taxrank = "Rank6")
#top.merge.Rank7 <- tax_glom(top.merge, taxrank = "Rank7")

# Export function
export_table <- function(phyloseq, saveTSV = F, returnTSV = F){
	table.all <- merge(t(otu_table(phyloseq)), tax_table(phyloseq), by = "row.names", all = TRUE)
	if(saveTSV){
		write.table(table.all, paste(deparse(substitute(phyloseq)), ".txt", sep = ""), sep="\t", row.names = F)
		table.all
	}
	if(returnTSV){
		return(table.all)
	}
}

export_table(top.merge.Rank2, saveTSV = T, returnTSV = F)
export_table(top.merge.Rank3, saveTSV = T, returnTSV = F)
export_table(top.merge.Rank4, saveTSV = T, returnTSV = F)
export_table(top.merge.Rank5, saveTSV = T, returnTSV = F)



