
# Setup -------------------------------------------------------------------

library(dplyr)
library(diversitree)

# Read in diversitree functions
source("code/diversitree_fns.R")
# Define function to round to specified multiple
mround <- function(x, base, FUN){
    base * FUN(x / base) 
}
# Define function to capitalise first letter
simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep = "", collapse = " ")
}

colourpal <- readRDS("data/misc/colourpal_mute.Rds")
birds_tree <- readRDS("data/raw/phy.Rds")
spp_colour <- 
    readRDS("data/sorted/summarise-colour/spp_colour.Rds") %>% 
    filter(
        patch == "crown",
        # Remove species with NA threat or trade status
        !(is.na(traded)),
        !(is.na(iucn)),
        !(is.na(CU_c_spMax))
    ) %>% 
    mutate(species = gsub("[.]", "_", species),
           traded = factor(traded,
                           levels = c("no", "yes"),
                           labels = c("not-traded", "traded")),
           threat = factor(iucn,
                           levels = c("LC", "NT", "VU", "EN", "CR"),
                           labels = c(rep("not-threatened",2), 
                                      rep("threatened",3))),
           threat_trade =
               factor(
                   case_when(traded == "traded" & threat == "threatened" ~ 
                                 "threatened_traded",
                             traded == "traded" & threat == "not-threatened" ~
                                 "not-threatened_traded",
                             traded == "not-traded" & threat == "threatened" ~
                                 "threatened_not-traded",
                             traded == "not-traded" & threat == "not-threatened" ~
                                 "not-threatened_not-traded"),
                   levels = c("threatened_traded",
                              "threatened_not-traded",
                              "not-threatened_traded",
                              "not-threatened_not-traded")))

# Create categories for uniqueness
unique_class <- classInt::classIntervals(spp_colour$CU_c_spMax, 4, style = "equal")

spp_colour <- 
    spp_colour %>%
    mutate(unique_class =
               cut(spp_colour$CU_c_spMax, 
                   unique_class$brks, 
                   include.lowest = TRUE))

# Get density values for histogram
unique_density <- density(spp_colour$CU_c_spMax)

# Get the hex value for each colour category
spp_colour$cat_hex <- 
    as.character(colourpal[match(spp_colour$colour_cat, names(colourpal))])

# Prune tree --------------------------------------------------------------

# Identify species in common between the tree and the colour spp_colour
sp_in_common <- 
    intersect(spp_colour$species, birds_tree$tip.label) 
# Prune the tree
birds_tree <- 
    drop.tip(birds_tree, 
             birds_tree$tip.label[-na.omit(match(sp_in_common, 
                                                 birds_tree$tip.label))])

# Sort --------------------------------------------------------------------
# Organize internal structure of phylogeny to get a ladderize effect
birds_tree <- ladderize(birds_tree, right = FALSE)
th <- max(branching.times(birds_tree))

# Get spp_colour rows in same order as they appear in the phylogeny
rownames(spp_colour) <- spp_colour$species
spp_colour <- spp_colour[birds_tree$tip.label,] 

# Define colour uniqueness colour palette, with grey for values less than 50
colour_unique_pal <- c("grey90", viridis::viridis(99))
# Cut colour uniqueness values
spp_colour_uniqueCuts <-
    cut(spp_colour$CU_c_spMax, 
        breaks = c(0, seq(50, max(spp_colour$CU_c_spMax), length.out = 100)))
# Name colour palette
names(colour_unique_pal) <- levels(spp_colour_uniqueCuts)
# Add to dataframe
spp_colour$unique_class_hist <- spp_colour_uniqueCuts
spp_colour$uniqueness_col <- colour_unique_pal[as.numeric(spp_colour_uniqueCuts)]
# Uniqueness from low (purple) to high (yellow)
birds_tree$tip.state <- spp_colour

# Also add to density df
unique_density_cuts <-
    cut(unique_density$x, breaks = c(0, seq(50, max(spp_colour$CU_c_spMax), length.out = 100)))
unique_density_cuts[which(is.na(unique_density_cuts))] <- levels(unique_density_cuts)[100]
# Match to colours
unique_density$col <- colour_unique_pal[unique_density_cuts]

# Order phylo ------------------------------------------------------------
ordered_tips <- 
    birds_tree$edge[birds_tree$edge[,2] <= length(birds_tree$tip.label), 2]

# Use rotateConstr to make sure tips are in order according to families
test <- cbind(spp_colour[ordered_tips,], seq = 1:length(ordered_tips))

fams <- unique(spp_colour$family)
famseq <- unique(test$family)

test2 <- test

for (i in 1:length(fams)) {
    stp1 <- test[which(test$family == famseq[i]),]
    # Force in sequence
    if (i == 1) {
        stp2 <- stp1
    }
    if (i > 1) {
        stp2 <- rbind(stp2, stp1)
    }
}
rownames(stp2) <- stp2$species

birds_tree <- rotateConstr(birds_tree, stp2$species)
birds_tree$tip.state <- stp2

## Force edges in sequence
t = 0
for (i in 1:length(birds_tree$edge[,2])) {
    if (birds_tree$edge[i,2] <= length(birds_tree$tip.label)) { # is tip?
        t <- t + 1
        birds_tree$edge[i,2] <- t
    }
    else{
        birds_tree$edge[i,2] <- birds_tree$edge[i,2]
    }
}

birds_tree$tip.label <- stp2$species

# Sequence black gray for families
plotord <- rev(unique(stp2$family))

if (length(plotord) == length(unique(spp_colour$family))) {   
    cat("Cool!!\nIt worked. \nfamily names match") 
} else {     
    cat("ops... something went wrong. Try again. \nReview family names...") 
}

# Write --------------------------------------------------------------------
xx <- c(rep(c("black", "grey"), floor(length(plotord)/2)))
# If length xx is even, add 'black'
if (length(xx) < length(plotord)) xx <- c(xx,'black')
xx <- xx[order(plotord)]

# Reduce palette to the colours actually present
colourpal <-
    colourpal[which(names(colourpal) %in% unique(spp_colour$colour_cat))]
# Sort legend labels
colour_labs <- 
    sapply(names(colourpal), simpleCap)
colour_labs <- gsub("_", " ", colour_labs)

# Save everything needed for analyses
saveRDS(birds_tree, "data/sorted/phylo/phylo.Rds")
