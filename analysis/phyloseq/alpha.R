# Alpha diversity over metadata variables 

#Colin Brislawn, 2015-12-22
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer", "scales")
library("vegan")
library("gridExtra")
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
# Like we said, all samples have good depth. We only lose one.

#####
# Go to folder for saving graphs
setwd("./alpha/")
list.files()



plot_richness_template <- function(phyloseq, xaxis, scalecolor, metrics, uselog=F, savegraph=T){
	plot <- plot_richness(phyloseq, x=xaxis, color=scalecolor, measures=metrics)
	plot$layers <- plot$layers[-1] # This pull off the geom_point layer, so we can add them back with different size.
	plot <- plot + geom_point(size=3, alpha=.3) + scale_colour_gradient(low = "lightblue", high = "navy", breaks = c(1, 0)) +
		stat_smooth(method=lm, color ="gray30", aes(fill = factor(PhoticZone))) + 
		scale_fill_manual(values = c("#222222", "khaki1"))
	#stat_smooth(method = "lm")
	
	if(uselog){
		plot <- plot + scale_x_log10()
		filename <- paste(xaxis, ".log10.png", sep = "")
	}
	
	if(savegraph){
		filename <- paste(xaxis, ".png", sep = "")
		ggsave(filename, width = 7, height = 4, units = 'in', dpi = 200, scale = 1.2)
	}
	return(plot)
}







#####
# Select one cohort
#cohort.rarefied <- subset_samples(filtered.rarefied)                      # All
#cohort.rarefied <- subset_samples(filtered, COHORT %in% c("12A", "12B"))  
cohort <- subset_samples(filtered, !is.na(MatDepth.cryosection.rel))       # Mat samples
cohort

# Add a factor for light and dark
tmp <- data.frame(sample_data(cohort))
tmp$PhoticZone <- ifelse(tmp$irradiance.PAR > 0, "Light", "Dark")
tmp$PhoticZone[is.na(tmp$PhoticZone)] <- "Dark"
head(tmp)
sample_data(cohort) <- tmp

# rarefy (subsample) (I'm sorry. This is needed when primarily using observed counts.) 
min(sample_sums(cohort))
cohort.rar <- rarefy_even_depth(cohort, sample.size = min(sample_sums(cohort)), replace = F, trimOTUs = F, rngseed = 711)
sort(sample_sums(cohort.rar))[1:5]

cohort.rar

# Posible measure are: c("Observed", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")




#Let’s build a 6 panel fig that shows just the observed OTUs plotted against 
#the following vars: Relative Depth (r.u.); %C Replaced (r.u.), Gross PS (µM O2 sec-1);
#Irradiance (µmol phot m-2 sec-1); dO2  (µM O2 ); Porosity (r.u.). 
#Can we lose the word “observed” over each plot? And… show the color map only once?
# Yes. Yes we can.
# Also, ask which ones of these should use log scale on x axis.
# Also, get R2 values and slopes for these regression overlays

metrics = c("Observed")
the_y_lim <- ylim(400,870) # For the 'Observed' metric


p1 <- plot_richness_template(cohort.rar, "MatDepth.cryosection.rel", "MatDepth.cryosection.rel", metrics, savegraph = F)
p1$layers <- p1$layers[1]
p1 <- p1 + theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none") +
	ylab("Observed OTUs") + xlab("Relative Depth (r.u.)") + stat_smooth() +
	annotate("text", .05, 850, label = "A", size = 6) 
p1

p2 <- plot_richness_template(cohort.rar, "pctC.replaced", "MatDepth.cryosection.rel", metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank()) +
	ylab("Observed OTUs") + xlab(bquote("%C Replaced ("*log[10]*" r.u.)")) + the_y_lim +
	annotate("text", 4, 850, label = "C", size = 6)
p2

p3 <- plot_richness_template(cohort.rar, "Gross.PS.rate", "MatDepth.cryosection.rel", metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank()) +
	scale_fill_manual(values = c("khaki1")) +
	ylab("Observed OTUs") + xlab(bquote('Gross PS (µM'*~O[2]~ s^-1*')')) + the_y_lim +
	annotate("text", 1, 850, label = "B", size = 6)
p3

p4 <- plot_richness_template(cohort.rar, "irradiance.PAR", "MatDepth.cryosection.rel", metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank(), plot.margin = unit(c(6,6,50,6), "pt")) +
	scale_fill_manual(values = c("khaki1")) +
	ylab("Observed OTUs") + xlab(bquote("Irradiance (µmol phot"*~m^-2 ~s^-1*")")) + the_y_lim +
	annotate("text", 5, 850, label = "E", size = 6)
p4

p5 <- plot_richness_template(cohort.rar, "dissolved.ox", "MatDepth.cryosection.rel", metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), 
				legend.position=c(-.8, -.45), legend.box = "horizontal", legend.direction = "horizontal", legend.background = element_blank(),
				axis.text.y=element_blank(), axis.title.y=element_blank(), plot.margin = unit(c(6,6,50,6), "pt")) +
	guides(colour = guide_colorbar("Relative Depth:"), fill = guide_legend("Photic Zone:")) +
	ylab("Observed OTUs") + xlab(bquote("d"*O[2]*" ("*log[10] ~µM ~O[2]*")")) +
	the_y_lim + scale_x_log10(breaks = c(250, 500, 1000, 1500)) +
	annotate("text", 1400, 850, label = "F", size = 6)
p5

p6 <- plot_richness_template(cohort.rar, "porosity", "MatDepth.cryosection.rel", metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none",
				plot.margin = unit(c(6,6,50,6), "pt")) +
	ylab("Observed OTUs") + xlab("Porosity (r.u.)") + the_y_lim +
	annotate("text", .2, 850, label = "D", size = 6)
p6

#g1
ggsave("6panel1.pdf", grid.arrange(p1, p3, p2, p6, p4, p5, ncol=3, widths=c(1.21,1,1), heights=c(1, 1.2)),
			 scale = 2.2, width = 89, height = 75, units = "mm", useDingbats=FALSE)
#


# Get R2 values for figure caption

r <- data.frame(sample_data(cohort.rar))
r$Observed <- as.data.frame(estimate_richness(cohort.rar, split = T, measures = "Observed"))
head(r)
names(r)
# Very long. I should put this in a loop or something.
# For each, I'll report R2 and slope

get_results <- function(thelm){
	test <- summary(thelm)
	output <- c(test$coefficients[2,1], #slope
	  test$coefficients[2,2], #Std. Error
	  test$coefficients[2,4], #p value
	  test$adj.r.squared
	)
	return(output)
}

get_results(lm(as.matrix(Observed) ~ Gross.PS.rate, r, PhoticZone=="Light"))

testresults <- data.frame(get_results(lm(as.matrix(Observed) ~ Gross.PS.rate, r, PhoticZone=="Light")))
testresults
names(testresults) <- "Gross.PS.rate.Light"
rownames(testresults) <- c("", "Slope", "Std. Error", "p value", "r squared")

testresults$log10.pctC.replaced.Light <- get_results(lm(as.matrix(Observed) ~ log10(pctC.replaced), r, PhoticZone=="Light"))
testresults$log10.pctC.replaced.Dark <- get_results(lm(as.matrix(Observed) ~ log10(pctC.replaced), r, PhoticZone=="Dark"))
testresults$porosity.Light <- get_results(lm(as.matrix(Observed) ~ (porosity), r, PhoticZone=="Light"))
testresults$porosity.Dark <- get_results(lm(as.matrix(Observed) ~ (porosity), r, PhoticZone=="Dark"))
testresults$irradiance.PAR.Light <- get_results(lm(as.matrix(Observed) ~ (irradiance.PAR), r, PhoticZone=="Light"))
testresults$log10.dissolved.ox.Light <- get_results(lm(as.matrix(Observed) ~ log10(dissolved.ox), r, PhoticZone=="Light"))
testresults$log10.dissolved.ox.Dark <- get_results(lm(as.matrix(Observed) ~ log10(dissolved.ox), r, PhoticZone=="Dark"))
testresults
write.table(t(testresults), "testresults.txt", F, F, sep = "\t")







# Another 6 panel figure with using the InvSimpson == Simpson's D

metrics = c("InvSimpson")
the_y_lim <- ylim(6,41) # For the 'InvSimpson' metric
the_y_lab <- ylab("Inverse Simpson")


p1 <- plot_richness_template(cohort.rar, "MatDepth.cryosection.rel", "MatDepth.cryosection.rel", metrics, savegraph = F)
p1$layers <- p1$layers[1]
p1 <- p1 + theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none") +
	the_y_lab + xlab("Relative Depth (r.u.)") + stat_smooth() + the_y_lim +
	annotate("text", .05, 38, label = "A", size = 6) 
p1

p2 <- plot_richness_template(cohort.rar, "pctC.replaced", "MatDepth.cryosection.rel", metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank()) +
	the_y_lab + xlab(bquote("%C Replaced ("*log[10]*" r.u.)")) + the_y_lim +
	annotate("text", 4, 38, label = "C", size = 6)
p2

p3 <- plot_richness_template(cohort.rar, "Gross.PS.rate", "MatDepth.cryosection.rel", metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank()) +
	scale_fill_manual(values = c("khaki1")) +
	the_y_lab + xlab(bquote('Gross PS (µM'*~O[2]~ s^-1*')')) + the_y_lim +
	annotate("text", 1, 38, label = "B", size = 6)
p3

p4 <- plot_richness_template(cohort.rar, "irradiance.PAR", "MatDepth.cryosection.rel", metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank(), plot.margin = unit(c(6,6,50,6), "pt")) +
	scale_fill_manual(values = c("khaki1")) +
	the_y_lab + xlab(bquote("Irradiance (µmol phot"*~m^-2 ~s^-1*")")) + the_y_lim +
	annotate("text", 5, 38, label = "E", size = 6)
p4

p5 <- plot_richness_template(cohort.rar, "dissolved.ox", "MatDepth.cryosection.rel", metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), 
				legend.position=c(-.8, -.45), legend.box = "horizontal", legend.direction = "horizontal", legend.background = element_blank(),
				axis.text.y=element_blank(), axis.title.y=element_blank(), plot.margin = unit(c(6,6,50,6), "pt")) +
	guides(colour = guide_colorbar("Relative Depth:"), fill = guide_legend("Photic Zone:")) +
	the_y_lab + xlab(bquote("d"*O[2]*" ("*log[10] ~µM ~O[2]*")")) +
	the_y_lim + scale_x_log10(breaks = c(250, 500, 1000, 1500)) +
	annotate("text", 1400, 38, label = "F", size = 6)
p5

p6 <- plot_richness_template(cohort.rar, "porosity", "MatDepth.cryosection.rel", metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none",
				plot.margin = unit(c(6,6,50,6), "pt")) +
	the_y_lab + xlab("Porosity (r.u.)") + the_y_lim +
	annotate("text", .2, 38, label = "D", size = 6)
p6

#g1
ggsave("6panel2.pdf", grid.arrange(p1, p3, p2, p6, p4, p5, ncol=3, widths=c(1.21,1,1), heights=c(1, 1.2)),
			 scale = 2.2, width = 89, height = 75, units = "mm", useDingbats=FALSE)
#









# 12 panel figure requested. Each of the 6 x-variables will be used, with
# 'richness' == observed and 'evenness' == Simpson's E stacked below.
# Panels A-L requested (6 x-axis * 2 metrics) but I think
# panels A-F (6 x-axis with two metrics) may be better. We'll see. 

# Set up special richness function which includes Simpsons Evenness.

estimate_richness_mod <- function(physeq, split=TRUE, measures=NULL){
	
	#' @importFrom vegan estimateR
	#' @importFrom vegan diversity
	#' @importFrom vegan fisher.alpha
	#' @export
	#' @examples 
	
	if( !any(otu_table(physeq)==1) ){
		# Check for singletons, and then warning if they are missing.
		# These metrics only really meaningful if singletons are included.
		warning(
			"The data you have provided does not have\n",
			"any singletons. This is highly suspicious. Results of richness\n",
			"estimates (for example) are probably unreliable, or wrong, if you have already\n",
			"trimmed low-abundance taxa from the data.\n",
			"\n",
			"We recommended that you find the un-trimmed data and retry."
		)
	}
	
	# If we are not splitting sample-wise, sum the species. Else, enforce orientation.
	if( !split ){
		OTU <- taxa_sums(physeq)		
	} else if( split ){
		OTU <- as(otu_table(physeq), "matrix")
		if( taxa_are_rows(physeq) ){ OTU <- t(OTU) }
	}
	
	# Define renaming vector:
	renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Pielou", "Simpson", "InvSimpson", "SimpsonD", "SimpsonE", "Fisher")
	names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", "pielou", "simpson", "invsimpson", "simpsond", "simpsone", "fisher")
	# If measures was not explicitly provided (is NULL), set to all supported methods
	if( is.null(measures) ){
		measures = as.character(renamevec)
	}
	# Rename measures if they are in the old-style
	if( any(measures %in% names(renamevec)) ){
		measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% measures]
	}
	
	# Stop with error if no measures are supported
	if( !any(measures %in% renamevec) ){
		stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
	}
	
	# Initialize to NULL
	outlist = vector("list")
	
	# Some standard diversity indices
	estimRmeas = c("Chao1", "Observed", "ACE")
	if( any(estimRmeas %in% measures) ){ 
		outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
	}
	if( "Shannon" %in% measures ){
		outlist <- c(outlist, list(shannon = diversity(OTU, index="shannon")))
	}
	if( "Pielou" %in% measures){
		#print("Starting Pielou")
		outlist <- c(outlist, list(pielou = diversity(OTU, index = "shannon")/log(specnumber(OTU))))
	}
	if( "Simpson" %in% measures ){
		outlist <- c(outlist, list(simpson = diversity(OTU, index="simpson")))
	}
	if( "InvSimpson" %in% measures ){
		outlist <- c(outlist, list(invsimpson = diversity(OTU, index="invsimpson")))
	}
	if( "SimpsonD" %in% measures ){
		# Different name but same metric as InvSimpson
		outlist <- c(outlist, list(simpsond = diversity(OTU, index="invsimpson")))
	}
	if( "SimpsonE" %in% measures ){
		#print("Starting SimpsonE")
		outlist <- c(outlist, list(simpsone = diversity(OTU, index="invsimpson")/specnumber(OTU)))
	}
	if( "Fisher" %in% measures ){
		fisher = tryCatch(fisher.alpha(OTU, se=TRUE), 
											warning=function(w){
												warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
												suppressWarnings(fisher.alpha(OTU, se=TRUE)[, c("alpha", "se")])
											}
		)
		if(!is.null(dim(fisher))){
			colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
			outlist <- c(outlist, list(fisher))
		} else {
			outlist <- c(outlist, Fisher=list(fisher))
		}
	}
	out = do.call("cbind", outlist)
	# Rename columns per renamevec
	namechange = intersect(colnames(out), names(renamevec))
	colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
	# Final prune to just those columns related to "measures". Use grep.
	colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out), ignore.case=TRUE)
	out = out[, sort(unique(unlist(colkeep))), drop=FALSE]
	# Make sure that you return a data.frame for reliable performance.
	out <- as.data.frame(out)
	return(out)
}


estimate_richness <- estimate_richness_mod

# This works:
head(estimate_richness_mod(cohort.rar))
# And so does this:
head(estimate_richness(cohort.rar))

# This should, but does not:
plot_richness(cohort.rar, x="MatDepth.cryosection.rel", color="MatDepth.cryosection.rel", measures=metrics)
# My guess is that it's using the estimate_richness from phyloseq, not the one
# I've made here. I can copy the full code into here and modify it outside of a
# library to make sure that it's using my estimate_richness function while, 
# graphing.
# Yep, that's what happening. I think the assignInNamespace() will let me hotwire this function. 
assignInNamespace("estimate_richness", estimate_richness_mod, ns = "phyloseq")
# God what a hack job. R has the structural integrity of a burning bridge of popsicle sticks. 
# OK, let's do this.

R.metrics = c("Observed")
E.metrics = c("SimpsonE")
R.the_y_lim <- ylim(400,870) # For the 'Observed' metric
R.the_y_lab <- ylab("Observed OTUs")
E.the_y_lim <- ylim(0.01, 0.05) # For the 'Simpson's E' metric
E.the_y_lab <- ylab("Evenness")


p.R.depth <- plot_richness_template(cohort.rar, "MatDepth.cryosection.rel", "MatDepth.cryosection.rel", R.metrics, savegraph = F)
p.R.depth$layers <- p.R.depth$layers[1]
p.R.depth <- p.R.depth + theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none",
															 axis.text.x=element_blank(), axis.title.x=element_blank()) +
	R.the_y_lab + xlab("Relative Depth (r.u.)") + stat_smooth() + R.the_y_lim
p.R.depth

p.E.depth <- plot_richness_template(cohort.rar, "MatDepth.cryosection.rel", "MatDepth.cryosection.rel", E.metrics, savegraph = F)
p.E.depth$layers <- p.E.depth$layers[1]
p.E.depth <- p.E.depth + theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none") +
	E.the_y_lab + xlab("Relative Depth (r.u.)") + stat_smooth() + E.the_y_lim
p.E.depth


p.R.ps <- plot_richness_template(cohort.rar, "Gross.PS.rate", "MatDepth.cryosection.rel", R.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank(),
				axis.text.x=element_blank(), axis.title.x=element_blank()) +
	scale_fill_manual(values = c("khaki1")) +
	R.the_y_lab + xlab(bquote('Gross PS (µM'*~O[2]~ s^-1*')')) + R.the_y_lim
p.R.ps

p.E.ps <- plot_richness_template(cohort.rar, "Gross.PS.rate", "MatDepth.cryosection.rel", E.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank()) +
	scale_fill_manual(values = c("khaki1")) +
	E.the_y_lab + xlab(bquote('Gross PS (µM'*~O[2]~ s^-1*')')) + E.the_y_lim
p.E.ps



p.R.replaced <- plot_richness_template(cohort.rar, "pctC.replaced", "MatDepth.cryosection.rel", R.metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank(),
				axis.text.x=element_blank(), axis.title.x=element_blank()) +
	R.the_y_lab + xlab(bquote("%C Replaced ("*log[10]*" r.u.)")) + R.the_y_lim
p.R.replaced

p.E.replaced <- plot_richness_template(cohort.rar, "pctC.replaced", "MatDepth.cryosection.rel", E.metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank()) +
	E.the_y_lab + xlab(bquote("%C Replaced ("*log[10]*" r.u.)")) + E.the_y_lim
p.E.replaced



p.R.porosity <- plot_richness_template(cohort.rar, "porosity", "MatDepth.cryosection.rel", R.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none",
				axis.text.x=element_blank(), axis.title.x=element_blank()) +
	R.the_y_lab + xlab("Porosity (r.u.)") + R.the_y_lim
p.R.porosity

p.E.porosity <- plot_richness_template(cohort.rar, "porosity", "MatDepth.cryosection.rel", E.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none",
				plot.margin = unit(c(6,6,50,6), "pt")) +
	E.the_y_lab + xlab("Porosity (r.u.)") + E.the_y_lim
p.E.porosity



p.R.irradiance <- plot_richness_template(cohort.rar, "irradiance.PAR", "MatDepth.cryosection.rel", R.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank(), 
				axis.text.x=element_blank(), axis.title.x=element_blank()) +
	scale_fill_manual(values = c("khaki1")) +
	R.the_y_lab + xlab(bquote("Irradiance (µmol phot"*~m^-2 ~s^-1*")")) + R.the_y_lim
p.R.irradiance

p.E.irradiance <- plot_richness_template(cohort.rar, "irradiance.PAR", "MatDepth.cryosection.rel", E.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank(),
				plot.margin = unit(c(6,6,50,6), "pt")) +
	scale_fill_manual(values = c("khaki1")) +
	E.the_y_lab + xlab(bquote("Irradiance (µmol phot"*~m^-2 ~s^-1*")")) + E.the_y_lim
p.E.irradiance



p.R.o2 <- plot_richness_template(cohort.rar, "dissolved.ox", "MatDepth.cryosection.rel", R.metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), 
				legend.position="none",
				axis.text.y=element_blank(), axis.title.y=element_blank(), 
				axis.text.x=element_blank(), axis.title.x=element_blank()) +
	guides(colour = guide_colorbar("Relative Depth:"), fill = guide_legend("Photic Zone:")) +
	R.the_y_lab + xlab(bquote("d"*O[2]*" ("*log[10] ~µM ~O[2]*")")) +
	R.the_y_lim + scale_x_log10(breaks = c(250, 500, 1000, 1500))
p.R.o2

p.E.o2 <- plot_richness_template(cohort.rar, "dissolved.ox", "MatDepth.cryosection.rel", E.metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), 
				legend.position=c(-.8, -.60), legend.box = "horizontal", legend.direction = "horizontal", legend.background = element_blank(),
				axis.text.y=element_blank(), axis.title.y=element_blank(), 
				plot.margin = unit(c(6,6,50,6), "pt"),
				panel.background = element_blank(),
				panel.ontop = F) +
	guides(colour = guide_colorbar("Relative Depth:"), fill = guide_legend("Photic Zone:")) +
	E.the_y_lab + xlab(bquote("d"*O[2]*" ("*log[10] ~µM ~O[2]*")")) +
	E.the_y_lim + scale_x_log10(breaks = c(250, 500, 1000, 1500))
p.E.o2


largepanel <- grid.arrange(p.R.depth, p.R.ps, p.R.replaced,
													 p.E.depth, p.E.ps, p.E.replaced, 
													 p.R.porosity, p.R.irradiance, p.R.o2,
													 p.E.porosity, p.E.irradiance, p.E.o2,
													 ncol=3, widths=c(1.22,1,1), heights=c(1, 1.4, 1, 1.8))
# I find that sometimes this fails to render with the error
# "polygon edge not found"
# A workaround is to change the width slightly (say by 0.01) and run grid.arrange again.
#large panel
ggsave("12panel thin.pdf", largepanel,
			 scale = 1.8, width = 89, height = 130, units = "mm", useDingbats=FALSE)
# 





# Rehashing of the 12 panel figure, for a horizontal layout


R.metrics = c("Observed")
E.metrics = c("SimpsonE")
R.the_y_lim <- ylim(400,870) # For the 'Observed' metric
R.the_y_lab <- ylab("Observed OTUs")
E.the_y_lim <- ylim(0.01, 0.05) # For the 'Simpson's E' metric
E.the_y_lab <- ylab("Evenness")


# The initial graphs are not different 
p.R.depth
p.E.depth

p.R.ps
p.E.ps

p.R.replaced
p.E.replaced


p.R.porosity <- plot_richness_template(cohort.rar, "porosity", "MatDepth.cryosection.rel", R.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none",
				axis.text.x=element_blank(), axis.title.x=element_blank(),
				axis.text.y=element_blank(), axis.title.y=element_blank()) +
	R.the_y_lab + xlab("Porosity (r.u.)") + R.the_y_lim
p.R.porosity

p.E.porosity <- plot_richness_template(cohort.rar, "porosity", "MatDepth.cryosection.rel", E.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none",
				axis.text.y=element_blank(), axis.title.y=element_blank()
				#plot.margin = unit(c(6,6,50,6), "pt")
				) +
	E.the_y_lab + xlab("Porosity (r.u.)") + E.the_y_lim
p.E.porosity



p.R.irradiance <- plot_richness_template(cohort.rar, "irradiance.PAR", "MatDepth.cryosection.rel", R.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank(), 
				axis.text.x=element_blank(), axis.title.x=element_blank()) +
	scale_fill_manual(values = c("khaki1")) +
	R.the_y_lab + xlab(bquote("Irradiance (µmol phot"*~m^-2 ~s^-1*")")) + R.the_y_lim
p.R.irradiance

p.E.irradiance <- plot_richness_template(cohort.rar, "irradiance.PAR", "MatDepth.cryosection.rel", E.metrics, savegraph = F) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), legend.position="none", 
				axis.text.y=element_blank(), axis.title.y=element_blank()
				#plot.margin = unit(c(6,6,50,6), "pt")
				) +
	scale_fill_manual(values = c("khaki1")) +
	E.the_y_lab + xlab(bquote("Irradiance (µmol phot"*~m^-2 ~s^-1*")")) + E.the_y_lim
p.E.irradiance



p.R.o2 <- plot_richness_template(cohort.rar, "dissolved.ox", "MatDepth.cryosection.rel", R.metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), 
				legend.position="right",
				axis.text.y=element_blank(), axis.title.y=element_blank(), 
				axis.text.x=element_blank(), axis.title.x=element_blank()) +
	guides(colour = guide_colorbar("Relative Depth:", reverse = T), fill = guide_legend("Photic Zone:")) +
	R.the_y_lab + xlab(bquote("d"*O[2]*" ("*log[10] ~µM ~O[2]*")")) +
	R.the_y_lim + scale_x_log10(breaks = c(250, 500, 1000, 1500))
p.R.o2

p.E.o2 <- plot_richness_template(cohort.rar, "dissolved.ox", "MatDepth.cryosection.rel", E.metrics, savegraph = F, uselog = T) +
	theme(strip.background = element_blank(),	strip.text.x = element_blank(), 
				legend.position= "right",
					#c(-.8, -.60), legend.box = "horizontal", legend.direction = "horizontal", legend.background = element_blank(),
				axis.text.y=element_blank(), axis.title.y=element_blank(), 
				#plot.margin = unit(c(6,6,50,6), "pt"),
				panel.background = element_blank(),
				panel.ontop = F) +
	guides(colour = guide_colorbar("Relative Depth:", reverse = T), fill = guide_legend("Photic Zone:")) +
	E.the_y_lab + xlab(bquote("d"*O[2]*" ("*log[10] ~µM ~O[2]*")")) +
	E.the_y_lim + scale_x_log10(breaks = c(250, 500, 1000, 1500))
p.E.o2


largepanel <- grid.arrange(p.R.depth, p.R.ps, p.R.replaced, p.R.porosity, p.R.irradiance, p.R.o2,
													 p.E.depth, p.E.ps, p.E.replaced, p.E.porosity, p.E.irradiance, p.E.o2,
													 ncol=6, widths=c(1.22,1,1, 1,1,1.8), heights=c(1, 1.2))
#largepanel
ggsave("12panel wide.pdf", largepanel,
			 scale = 2, width = 183, height = 60, units = "mm", useDingbats=FALSE)
# 





######
# 
r <- data.frame(sample_data(cohort.rar))
r$Observed <- as.data.frame(estimate_richness_mod(cohort.rar, split = T, measures = "Observed"))
r$SimpsonE <- as.data.frame(estimate_richness_mod(cohort.rar, split = T, measures = "SimpsonE"))
head(r)
names(r)
# Let's make a loop, just like above
# For each, I'll report R2 and slope

get_results <- function(thelm){
	test <- summary(thelm)
	output <- c(test$coefficients[2,1], #slope
							test$coefficients[2,2], #Std. Error
							test$coefficients[2,4], #p value
							test$adj.r.squared
	)
	return(output)
}

summary(lm(as.matrix(Observed) ~ Gross.PS.rate, r, PhoticZone=="Light"))
get_results(lm(as.matrix(Observed) ~ Gross.PS.rate, r, PhoticZone=="Light"))

testresults <- data.frame(get_results(lm(as.matrix(Observed) ~ Gross.PS.rate, r, PhoticZone=="Light")))
testresults
names(testresults) <- "observed.gross.PS.rate.Light"
rownames(testresults) <- c("Slope", "Std. Error", "p value", "r squared")

# Observed (richness)
testresults$observed.log10.pctC.replaced.Light <- get_results(lm(as.matrix(Observed) ~ log10(pctC.replaced), r, PhoticZone=="Light"))
testresults$observed.log10.pctC.replaced.Dark <- get_results(lm(as.matrix(Observed) ~ log10(pctC.replaced), r, PhoticZone=="Dark"))
testresults$observed.porosity.Light <- get_results(lm(as.matrix(Observed) ~ (porosity), r, PhoticZone=="Light"))
testresults$observed.porosity.Dark <- get_results(lm(as.matrix(Observed) ~ (porosity), r, PhoticZone=="Dark"))
testresults$observed.irradiance.PAR.Light <- get_results(lm(as.matrix(Observed) ~ (irradiance.PAR), r, PhoticZone=="Light"))
testresults$observed.log10.dissolved.ox.Light <- get_results(lm(as.matrix(Observed) ~ log10(dissolved.ox), r, PhoticZone=="Light"))
testresults$observed.log10.dissolved.ox.Dark <- get_results(lm(as.matrix(Observed) ~ log10(dissolved.ox), r, PhoticZone=="Dark"))

# SimpsonE (evenness)
get_results(lm(as.matrix(SimpsonE) ~ Gross.PS.rate, r, PhoticZone=="Light"))

testresults$evenness.Gross.PS.rate.Light <- get_results(lm(as.matrix(SimpsonE) ~ Gross.PS.rate, r, PhoticZone=="Light"))
testresults$evenness.log10.pctC.replaced.Light <- get_results(lm(as.matrix(SimpsonE) ~ log10(pctC.replaced), r, PhoticZone=="Light"))
testresults$evenness.log10.pctC.replaced.Dark <- get_results(lm(as.matrix(SimpsonE) ~ log10(pctC.replaced), r, PhoticZone=="Dark"))
testresults$evenness.porosity.Light <- get_results(lm(as.matrix(SimpsonE) ~ (porosity), r, PhoticZone=="Light"))
testresults$evenness.porosity.Dark <- get_results(lm(as.matrix(SimpsonE) ~ (porosity), r, PhoticZone=="Dark"))
testresults$evenness.irradiance.PAR.Light <- get_results(lm(as.matrix(SimpsonE) ~ (irradiance.PAR), r, PhoticZone=="Light"))
testresults$evenness.log10.dissolved.ox.Light <- get_results(lm(as.matrix(SimpsonE) ~ log10(dissolved.ox), r, PhoticZone=="Light"))
testresults$evenness.log10.dissolved.ox.Dark <- get_results(lm(as.matrix(SimpsonE) ~ log10(dissolved.ox), r, PhoticZone=="Dark"))

# looking into negative adjusted R2 value
get_results(lm(as.matrix(SimpsonE) ~ log10(pctC.replaced), r, PhoticZone=="Dark"))
summary(lm(as.matrix(SimpsonE) ~ log10(pctC.replaced), r, PhoticZone=="Dark"))

testresults
dim(testresults)
write.table(t(testresults), "testresults.txt", F, F, sep = "\t")





