# Analyze community composition of NEON and DoB soil fungi

RAREFACTION_ITS <- 10000
MIN_SEQ_DEPTH <- 2500

library(phyloseq)

neon <- readRDS("./code/NEON_ITS_phyloseq_DL08-13-2019.Rds")
nrow(otu_table(neon))
sum(sample_sums(neon)<MIN_SEQ_DEPTH) # This many samples will be removed at MIN_SEQ_DEPTH

hist(sample_sums(neon), breaks=60, main=paste("Sequencing depth of NEON ITS samples,\nn =", length(sample_sums(neon))))
hist(sample_sums(neon), breaks=60, main=paste("Sequencing depth of NEON ITS samples\nwith >", MIN_SEQ_DEPTH, "reads, n =", nrow(otu_table(neon)) - sum(sample_sums(neon)<MIN_SEQ_DEPTH)))

# # rarefaction
# # actually, don't rarefy for composition (following recommendations by 
# # McMurdie & Holmes, 2014)
# sum(sample_sums(neon) < RAREFACTION_ITS) # this many samples would be lost at current rarefaction level
# neon_rare <- rarefy_even_depth(physeq=neon,sample.size=RAREFACTION_ITS,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)

# Aggregate to genus or family
neon_genus <- tax_glom(neon, taxrank = "Genus")
neon_family <- tax_glom(neon, taxrank = "Family")
# Note: this removes sequences with 'NA' genus or family,
# so it is possible to end up with empty samples. Remove these
# (and optionally other low-count samples):
which(sample_sums(neon_genus)==0)
which(sample_sums(neon_family)==0)
neon_genus <- prune_samples(sample_sums(neon_genus)>0, neon_genus)
neon_family <- prune_samples(sample_sums(neon_family)>0, neon_family)

neon_genus_mindepth <- prune_samples(sample_sums(neon_genus)>MIN_SEQ_DEPTH, neon_genus)

neon_mindepth <- prune_samples(sample_sums(neon)>MIN_SEQ_DEPTH, neon)

# Remove taxa with trivially large C.V., apply filter to keep only
# taxa which show up at least twice in at least two different samples.
# (Removes ~81% of all taxa in non-aggregated dataset
neon_mindepth_filttaxa <- filter_taxa(neon_mindepth, function(x) { sum(x >= 2) >= 2 }, TRUE)
# How much did this affect distribution of C.V.?
cv_nofilt <- apply(otu_table(neon_mindepth), 2, function(x) { sd(x) / mean(x) }) # before
cv_filttaxa <- apply(otu_table(neon_mindepth_filttaxa), 2, function(x) { sd(x) / mean(x) }) # after
hist(cv_nofilt) # before
hist(cv_filttaxa) # after

# Normalize counts by proportion of each sample's library size
neon_genus_prop <- transform_sample_counts(neon_genus, function(OTU) OTU/sum(OTU))
neon_family_prop <- transform_sample_counts(neon_family, function(OTU) OTU/sum(OTU))

neon_genus_mindepth_prop <- transform_sample_counts(neon_genus_mindepth, function(OTU) OTU/sum(OTU))

neon_mindepth_prop <- transform_sample_counts(neon_mindepth, function(OTU) OTU/sum(OTU))
neon_mindepth_filttaxa_prop <- transform_sample_counts(neon_mindepth_filttaxa, function(OTU) OTU/sum(OTU))

# Ordinate
neon_genus_ord <- ordinate(neon_genus_prop, "NMDS", "bray")
# p1 = plot_ordination(neon_genus, neon_genus_ord, type="taxa", color="Phylum", title="taxa")
p1 = plot_ordination(neon_genus_prop, neon_genus_ord, type="samples", title="Samples, aggregated to genus")
print(p1)

neon_family_ord <- ordinate(neon_family_prop, "NMDS", "bray")
p2 = plot_ordination(neon_family, neon_family_ord, type="samples", title="Samples, aggregated to family")
print(p2)

neon_genus_mindepth_ord <- ordinate(neon_genus_mindepth_prop, "NMDS", "bray")
p3 = plot_ordination(neon_genus_mindepth_prop, neon_genus_mindepth_ord, type="samples", 
                     color="domainID", title=paste0("Samples, aggregated to genus, mindepth=",MIN_SEQ_DEPTH))
print(p3)

p3.1 = plot_ordination(neon_genus_mindepth_prop, neon_genus_mindepth_ord, type="samples", 
                       color="soilMoisture", title=paste0("Samples, aggregated to genus, mindepth=",MIN_SEQ_DEPTH))
print(p3.1)

neon_mindepth_ord <- ordinate(neon_mindepth_prop, "NMDS", "bray")
p4 = plot_ordination(neon_mindepth_prop, neon_mindepth_ord, type="samples", 
                       color="domainID", title=paste0("Samples, mindepth=",MIN_SEQ_DEPTH))
print(p4)

neon_mindepth_filttaxa_ord <- ordinate(neon_mindepth_filttaxa_prop, "NMDS", "bray")
p5 <- plot_ordination(neon_mindepth_filttaxa_prop, neon_mindepth_filttaxa_ord, type="samples", 
                      color="domainID", title=paste0("Samples, filtered taxa, mindepth=",MIN_SEQ_DEPTH))
print(p5)
p5_outliers <- rownames(scores(neon_mindepth_filttaxa_ord)[which(scores(neon_mindepth_filttaxa_ord)[,1] < -0.4),])

# Most abundant taxa in outlier samples (normalized by number of samples)
head(sort(taxa_sums(prune_samples(p5_outliers, neon_mindepth_filttaxa_prop)), decreasing = TRUE))/length(p5_outliers)
# Most abundant taxa in non-outlier samples (normalized by number of samples)
head(sort(taxa_sums(prune_samples(!(rownames(otu_table(neon_mindepth_filttaxa_prop)) %in% p5_outliers), neon_mindepth_filttaxa_prop)), decreasing = TRUE))/(length(sample_sums(neon_mindepth_filttaxa_prop))-length(p5_outliers))
# Within outliers, what is the abundance of the taxa that are 
# most abundant in non-outliers?
head(taxa_sums(prune_samples(p5_outliers, neon_mindepth_filttaxa_prop)), n=200)/length(p5_outliers)
# Answer: not abundant at all...

# Let's calculate the following:
# How many samples n have a sample-sum of the first k most
# abundant taxa that is equal to zero?
samples_with_zero_of_first_k_species <- function(physeq, k) {
  presence_absence <- as(object = otu_table(physeq), Class = "matrix") != 0
  # initial_boolean_vector <- rep(TRUE, nrow(presence_absence))
  cumulative_absence <- matrix(nrow=nrow(presence_absence), ncol=k, 
                               dimnames=list(rownames(presence_absence), colnames(presence_absence)[1:k]))
  cumulative_absence[,1] <- presence_absence[,1]
  for(i in 2:k) {
    cumulative_absence[,i] <- cumulative_absence[,i-1] | presence_absence[,i]
  }
  return(list(
    n = sum(!cumulative_absence[,k]),
    n_cumulative = apply(!cumulative_absence, 2, sum),
    sample_names = rownames(cumulative_absence)[which(cumulative_absence[,k]==FALSE)])
  )
}

samples_with_zero_of_first_k_species(neon_mindepth_filttaxa, 1499)
p5_outliers
# THESE MATCH EXACTLY
# This means that for some reason, the 12 outlier samples are outliers
# because they don't contain ANY trace of the top 1499 most abundant ASVs
# in the dataset.
# In fact, they're such outliers that we only had to set k = 190 to wittle
# down the sample names to this group, and then this group didn't dissolve 
# until k = 1500.
# I don't know why this was the case, but I believe this justifies removal
# of these outliers.

neon_mindepth_filttaxa_prop_rmoutliers <- prune_samples(!(rownames(otu_table(neon_mindepth_filttaxa_prop)) %in% p5_outliers), neon_mindepth_filttaxa_prop)

p6_ord <- ordinate(neon_mindepth_filttaxa_prop_rmoutliers, "NMDS", "bray")
p6 <- plot_ordination(neon_mindepth_filttaxa_prop_rmoutliers, p6_ord, type="samples", 
                      color="domainID", title=paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust = 0.5))
print(p6)

p6.1 <- plot_ordination(neon_mindepth_filttaxa_prop_rmoutliers, p6_ord, type="samples", 
                        color="horizon", title=paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust = 0.5))
p6.1

p6.2 <- plot_ordination(neon_mindepth_filttaxa_prop_rmoutliers, p6_ord, type="samples", 
                        color="soilInWaterpH", title=paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust = 0.5))
p6.2

p6.3 <- plot_ordination(neon_mindepth_filttaxa_prop_rmoutliers, p6_ord, type="samples", 
                        color="decimalLatitude", title=paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust = 0.5))
p6.3

p6.4 <- plot_ordination(neon_mindepth_filttaxa_prop_rmoutliers, p6_ord, type="samples", 
                        color="soilMoisture", title=paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust = 0.5))
p6.4

p6.5 <- plot_ordination(neon_mindepth_filttaxa_prop_rmoutliers, p6_ord, type="samples", 
                        color="elevation", title=paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust = 0.5))
p6.5

# This one doesn't work for some reason. Date is interpreted as character.
p6.6 <- plot_ordination(neon_mindepth_filttaxa_prop_rmoutliers, p6_ord, type="samples", 
                        color="collectDate.x", title=paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust = 0.5))
p6.6

neon_rmoutliers_scores <- as.data.frame(scores(p6_ord))
neon_rmoutliers_scores$sample <- rownames(neon_rmoutliers_scores)
neon_rmoutliers_sample_data <- as(sample_data(neon_mindepth_filttaxa_prop_rmoutliers), Class="data.frame")
neon_rmoutliers_sample_data$sample <- rownames(neon_rmoutliers_sample_data)
neon_rmoutliers_sample_data$collectDate <- as.Date(neon_rmoutliers_sample_data$collectDate.x)

neon_rmoutliers_sample_data %>%
  full_join(neon_rmoutliers_scores, by="sample") ->
  neon_rmoutliers_sample_data

ggplot(neon_rmoutliers_sample_data,
       aes(x=NMDS1, y=NMDS2, col=collectDate)) +
  geom_point() +
  ggtitle(paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust=0.5))
# no effect of collection date on composition at macro spatial scale






################
# Get climate data for NEON field sites

library(dplyr)
library(tidyr)
neon_sites0 <- read.csv("field-sites.csv", header=TRUE)
names(neon_sites0) <- tolower(names(neon_sites0))
neon_sites <- tidyr::extract(neon_sites0, "lat..long.", c("lat","lon"), "(-?[0-9]+.[0-9]+),\\s*(-?[0-9]+.[0-9]+)")
neon_sites %>%
  mutate(lat = as.numeric(lat),
         lon = as.numeric(lon)) %>%
  
  # Also filter for terrestrial sites only
  filter(site.type %in% c("Core Terrestrial", "Relocatable Terrestrial")) ->
  neon_sites

library(raster)
library(sp)
# Load current (actually 1970-2000 average) climate data
r_current <- getData("worldclim",var="bio",res=10)
r_current <- r_current[[c(1,12)]]
names(r_current) <- c("temp0","prec0")

# Get coordinates for NEON sites
lats <- neon_sites$lat
lons <- neon_sites$lon
coords <- data.frame(x=lons,y=lats)
points <- SpatialPoints(coords, proj4string = r_current@crs)
values_current <- extract(r_current, points)

cbind.data.frame(
  coordinates(points),
  values_current,
  siteID = neon_sites$site.id
) %>%
  mutate(temp0 = temp0/10,     # For some reason, Wordclim data comes
         prec0 = prec0/10) %>% #   with values multiplied by 10)
  rename(mat_celsius = temp0,
         map_cm = prec0) ->
  clim_df

# Now remap with MAT and MAP
neon_rmoutliers_sample_data %>%
  left_join(clim_df, by="siteID") ->
  neon_rmoutliers_sample_data
  
# MAT
ggplot(neon_rmoutliers_sample_data,
       aes(x=NMDS1, y=NMDS2, col=mat_celsius)) +
  geom_point() +
  ggtitle(paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_colour_gradient(low = "#0000FF", high = "#FF0000",
                        na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour")

# MAP
ggplot(neon_rmoutliers_sample_data,
       aes(x=NMDS1, y=NMDS2, col=map_cm)) +
  geom_point() +
  ggtitle(paste0("Samples, mindepth=",MIN_SEQ_DEPTH,", filtered taxa,\n removed outliers")) +
  theme(plot.title = element_text(hjust=0.5))
