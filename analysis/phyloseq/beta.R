# Beta diversity, and finally, biplots

#Colin Brislawn, 2016-02-04

library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer", "scales")
library("vegan")
theme_set(theme_bw())
library("plyr")
library("reshape2")

#####
# Importing data

no_meta <- import_biom(file.path('../../data/otus_vsearch/otu_table_w_tax.biom'),
											 file.path('../../data/otus_vsearch/rep_set.tre'))
no_meta

meta <- import_qiime_sample_data(file.path('../../metadata/metadata.txt'))

full <- merge_phyloseq(meta, no_meta)
full



#####

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


############
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
# All samples have good depth. We only lose one.

# rarefy (subsample)
min(sample_sums(filtered))
filtered.rarefied <- rarefy_even_depth(filtered, replace = F, trimOTUs = F, rngseed = 711)
sort(sample_sums(filtered.rarefied))[1:5]

filtered.rarefied

# Select right cohort (no mock communities)
cohort.rarefied <- subset_samples(filtered.rarefied, !is.na(MatDepth.cryosection.rel))       # Mat samples
cohort.rarefied
sort(sample_sums(cohort.rarefied))[1:5]

# Convert NA to zeros in the sample_data
sample_data(cohort.rarefied)[is.na(sample_data(cohort.rarefied))] <- 0
sample_data(cohort.rarefied)

#####
# Go to folder for saving graphs
setwd("./beta/")
list.files()




#####
# Initial Beta diversity plots on merged data

# Use the filtered, rarefied cohort
cohort.rarefied

cohort.bray = ordinate(cohort.rarefied, method = "PCoA", distance = "bray")
cohort.unifrac = ordinate(cohort.rarefied, method = "PCoA", "unifrac")
cohort.wunifrac = ordinate(cohort.rarefied, method = "PCoA", "unifrac", weighted = T)

# Option to save precalculated distances. 
#save.image("beta.RData.gz", compress = T)
#load("beta.RData.gz")

cohort.bray.plot <- plot_ordination(cohort.rarefied, cohort.bray, 
																		color = "MatDepth.cryosection.rel",
																		#shape = "Experiment",
																		title = "PCoA Bray-Curtis")
cohort.bray.plot

cohort.unifrac.plot <- plot_ordination(cohort.rarefied, cohort.unifrac, 
																			 color = "MatDepth.cryosection.rel", 
																			 #shape = "Experiment",
																			 title = "PCoA unweighted UniFrac")
cohort.unifrac.plot 

cohort.wunifrac.plot <- plot_ordination(cohort.rarefied, cohort.wunifrac, 
																				color = "MatDepth.cryosection.rel",
																				#shape = "Tree",
																				title = "PCoA weighted UniFrac")
cohort.wunifrac.plot

#Weighted and unweighted are both interesting. Weighted UniFrac will be the metric of choice.

ggsave('beta bray.png', cohort.bray.plot, width = 5, height = 3, units = 'in', dpi = 200, scale = 1.4)
ggsave('beta unifrac.png', cohort.unifrac.plot, width = 5, height = 3, units = 'in', dpi = 200, scale = 1.4)
ggsave('beta wunifrac.png', cohort.wunifrac.plot, width = 5, height = 3, units = 'in', dpi = 200, scale = 1.4)



#####
# adding taxa to metadata for resolution in biplot
# Add the most common taxa (at a given level) to the sample_data() so that
# these taxa can be graphed as be vectors.

glom5 <- tax_glom(cohort.rarefied, taxrank = "Rank5")
#glom4 <- tax_glom(glom5, taxrank = "Rank4")
glom <- glom5
head(tax_table(glom))

glom.top <- prune_taxa(names(sort(taxa_sums(glom), TRUE))[0:10], glom)
glom.top
sum(taxa_sums(glom.top)) / sum(taxa_sums(glom))

tax_table(glom.top)
# These are the same OTUs which Hans sent me. Good.

glom.top.table <- otu_table(glom.top)
glom.top.table
taxa_sums(glom.top.table)
row.names(glom.top.table)
tax_table(glom.top)

sort(taxa_sums(glom.top.table), T)

# Manual annotation of row names. 
rownames(glom.top.table) <- c("o__Spirochaetales f__Spirochaetaceae", 
															"o__Oscillatoriales f__Phormidiaceae",
															"o__Pseudanabaenales f__Pseudanabaenaceae",
															"o__SBR1031 f__A4b",
															"o__Rhodospirillales",
															"o__Rhodobacterales f__Rhodobacteraceae",
															"o__Rhizobiales",
															"o__Rhodospirillales f__Rhodospirillaceae",
															"o__Bacteroidales f__ML635J-40",
															"o__Phycisphaerales")
head(t(glom.top.table))
head(sample_data(cohort.rarefied))


cohort.rarefied.w_tax <- cohort.rarefied

tmp <- merge(sample_data(cohort.rarefied), t(glom.top.table), by = "row.names")
head(tmp)
row.names(tmp) <- tmp$Row.names

# Quick graph of tmp to show off species abundance
ggplot(tmp, aes(MatDepth.cryosection.rel, o__Phycisphaerales)) + geom_point()
tmp.short <- tmp
tmp.short$BarcodeSequence <- NULL
tmp.short$LinkerPrimerSequence <- NULL
tmp.short$InputFileName <- NULL
tmp.short$HansSampleID <- NULL
tmp.short$MatDepth.cryosection.mm <- NULL
tmp.short$MatDepth.cryosection.mm.scaled <- NULL
tmp.short$Cryosection.replicate <- NULL
tmp.short$Gross.PS.rate <- NULL
tmp.short$pctC.replaced <- NULL
tmp.short$pctC.replaced.perhr <- NULL
tmp.short$X13C.isotoperatio.rate <- NULL
tmp.short$irradiance.PAR <- NULL
tmp.short$dissolved.ox <- NULL
tmp.short$porosity <- NULL
tmp.short$Description <- NULL
tmp.short$Row.names <- NULL
tmp.short$X.SampleID <- NULL
# Also remove uninteresting taxa for this graph.
# Actually, keep these for new graph.
# tmp.short$`o__Spirochaetales f__Spirochaetaceae` <- NULL
# tmp.short$`o__SBR1031 f__A4b` <- NULL
# tmp.short$`o__Rhodospirillales f__Rhodospirillaceae` <- NULL
# tmp.short$`o__Rhodobacterales f__Rhodobacteraceae` <- NULL
# tmp.short$o__Rhizobiales <- NULL
# tmp.short$o__Rhodospirillales <- NULL
# tmp.short$o__Phycisphaerales <- NULL

head(tmp.short)

tmp.melt <- melt(tmp.short, id=c("MatDepth.cryosection.rel"))
head(tmp.melt)
table(tmp.melt$variable)

g <- ggplot(tmp.melt, aes(MatDepth.cryosection.rel, value, color=MatDepth.cryosection.rel))
g <- g + geom_point() + facet_wrap(~variable, ncol = 2) + 
	scale_y_log10() +
	scale_colour_gradient(low = "lightblue", high = "navy", guide = guide_colorbar(title = "Relative Mat Depth", reverse = T)) +
	ylab("Rarified Abundance at the Family Level") +
	xlab("Relative Mat Depth")
g	
#ggsave("select families.png", g, scale = 1.5, units = "in", width = 5, height = 4)
ggsave("wunifrac CAP families.png", g, scale = 1.2, units = "mm", width = 183, height = 169)
ggsave("wunifrac CAP families.pdf", g, scale = 1.2, units = "mm", width = 183, height = 169, useDingbats = F)
sample_data(cohort.rarefied.w_tax) <- tmp
head(sample_data(cohort.rarefied.w_tax))
#







#####
# Biplots!

# capscale example
data(varespec)
head(varespec)
class(varespec)

data(varechem)
head(varechem)
class(varechem)
## Basic Analysis
vare.cap <- capscale(varespec ~ N + P + K + Condition(Al), varechem,
										 dist="bray")
vare.cap
plot(vare.cap)
#anova(vare.cap)



# similar thing using capscale.phyloseq
cohort.distance <- phyloseq::distance(cohort.rarefied.w_tax, method = "unifrac", type = "samples", weighted = T)
hist(cohort.distance)
sample_variables(cohort.rarefied.w_tax)
# Plot CAP with capscale 
vare.cap <- capscale(cohort.distance ~ irradiance.PAR + Gross.PS.rate + pctC.replaced + dissolved.ox + porosity +
										 	o__Spirochaetales.f__Spirochaetaceae +
										 	o__Oscillatoriales.f__Phormidiaceae +
										 	o__Pseudanabaenales.f__Pseudanabaenaceae +
										 	o__SBR1031.f__A4b +
										 	o__Rhodospirillales +
										 	o__Rhodobacterales.f__Rhodobacteraceae +
										 	o__Rhizobiales +
										 	o__Rhodospirillales.f__Rhodospirillaceae +
										 	o__Bacteroidales.f__ML635J.40 +
										 	o__Phycisphaerales,
										 data.frame(sample_data(cohort.rarefied.w_tax)))
plot(vare.cap)
#


# People like circles around the db-RDA vectors. So, here, have a circle:
# From http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2 
circleFun <- function(center = c(0,0), radius = 1, npoints = 100){
	r = radius
	tt <- seq(0, 2*pi, length.out = npoints)
	xx <- center[1] + r * cos(tt)
	yy <- center[2] + r * sin(tt)
	return(data.frame(x = xx, y = yy))
}

#geom_path will do open circles, geom_polygon will do filled circles
#plot + geom_path(aes(x = x, y=y), data = circleFun(c(0,0), 1, 100), color="grey")




######
# CAP / db-RDA
######
cohort.distance.ord = ordinate(cohort.rarefied.w_tax, "CAP", cohort.distance, ~ MatDepth.cryosection.rel + irradiance.PAR + Gross.PS.rate + pctC.replaced + dissolved.ox + porosity + 
															 	o__Spirochaetales.f__Spirochaetaceae +
															 	o__Oscillatoriales.f__Phormidiaceae +
															 	o__Pseudanabaenales.f__Pseudanabaenaceae +
															 	o__SBR1031.f__A4b +
															 	o__Rhodospirillales +
															 	o__Rhodobacterales.f__Rhodobacteraceae +
															 	o__Rhizobiales +
															 	o__Rhodospirillales.f__Rhodospirillaceae +
															 	o__Bacteroidales.f__ML635J.40 +
															 	o__Phycisphaerales
															 )
#plot_ordination(cohort.rarefied, cohort.distance.ord, "samples", color="MatDepth.cryosection.rel")
#plot(cohort.distance.ord)

# pull out elements of the biplot inside the ordination object
names(cohort.distance.ord)
cohort.distance.ord$CCA$biplot
vectors <- data.frame(cohort.distance.ord$CCA$biplot)
# The find the length (magnitude) of the vectors in the first two dimensions
vectors$length1and2 <-sqrt(vectors$CAP1^2+vectors$CAP2^2)
vectors

# Better names
row.names(data.frame(cohort.distance.ord$CCA$biplot))
vectornames <- c("Relative Depth", "Irradiance", "Gross PS", "%C Replaced", "Dissolved Oxygen", "Porosity",
								 #"o__Spirochaetales \nf__Spirochaetaceae", 
								 "f__Spirochaetaceae", 
								 #"o__Oscillatoriales \nf__Phormidiaceae",
								 "f__Phormidiaceae",
								 #"o__Pseudanabaenales \nf__Pseudanabaenaceae",
								 "f__Pseudanabaenaceae",
								 #"o__SBR1031 \nf__A4b",
								 "SBR1031 A4b",
								 "o__Rhodospirillales",
								 #"o__Rhodobacterales \nf__Rhodobacteraceae",
								 "f__Rhodobacteraceae",
								 "o__Rhizobiales",
								 #"o__Rhodospirillales \nf__Rhodospirillaceae",
								 "f__Rhodospirillaceae",
								 #"o__Bacteroidales \nf__ML635J-40",
								 "Bacteroidales ML635J-40",
								 "o__Phycisphaerales")

# Add vector categories 
vector.data <- data.frame(cohort.distance.ord$CCA$biplot)
vector.data$category <- c("Environment", "Environment",
													"Productivity", "Productivity",
													"Environment","Environment", 
													"Taxa", "Taxa", "Taxa", "Taxa", "Taxa",
													"Taxa", "Taxa", "Taxa", "Taxa", "Taxa")
envcolor = "#777777"
procolor = "#FF523F"
taxcolor = "#5CA9E0"
# For pallet creation: https://coolors.co/app/ff523f-777777-5ca9e0-1e3231-485665 
vector.data$color <- c(envcolor, envcolor,
											 procolor, procolor,
											 envcolor, envcolor, 
											 taxcolor, taxcolor, taxcolor, taxcolor, taxcolor,
											 taxcolor, taxcolor, taxcolor, taxcolor, taxcolor)

# Now Let's get those vectors back into the ggplot
plot <- plot_ordination(cohort.rarefied, cohort.distance.ord, "samples", color="MatDepth.cryosection.rel")
#plot
plot <- plot + 
	geom_path(aes(x = x, y=y), data = circleFun(c(0,0), 1, 100), color="lightgrey") + 
	geom_segment(data = vector.data, aes(x = 0, xend = CAP1, y = 0, yend = CAP2, color = category),
						 arrow = arrow(length = unit(0.25, "cm")), color = vector.data$color) + 
	geom_text(data=data.frame(cohort.distance.ord$CCA$biplot), 
						aes(x=CAP1, y=CAP2, label=vectornames), 
						size = 4, 
						color = vector.data$color
						#color="black"
						) +
	scale_colour_gradient(low = "lightblue", high = "navy", breaks = c(1, 0),
												guide = guide_colorbar(title = "Relative Depth", reverse = T)) +
	theme(legend.position = c(.85, .2)) 
plot
#ggsave("wunifrac CAP biplot tax.png", plot, scale = 1.2, width = 183, height = 183, units = 'mm', dpi = 200)
ggsave("wunifrac CAP biplot.pdf", plot, scale = 1, width = 183, height = 183, units = 'mm', dpi = 200, useDingbats=FALSE)
ggsave("wunifrac CAP biplot bigger.pdf", plot, scale = 1.6, width = 89, height = 86, units = 'mm', dpi = 200, useDingbats=FALSE)
#
# two column figures are 178 mm
# max depth is 210 mm





######
# Additional capscale / CAP / db-RDA examples
## Avoid negative eigenvalues with additive constant
capscale(varespec ~ N + P + K + Condition(Al), varechem, dist="bray", add =TRUE)
## Avoid negative eigenvalues by taking square roots of dissimilarities
capscale(varespec ~ N + P + K + Condition(Al), varechem,
				 dist = "bray", sqrt.dist= TRUE)
## Principal coordinates analysis with extended dissimilarities
capscale(varespec ~ 1, dist="bray", metaMDS = TRUE)




