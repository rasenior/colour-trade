
# Setup ------------------------------------------------------------------------

library(dplyr)
library(geometry)

# Helper functions --------------------------------------------------------

mround <- function(x,base,FUN){
    base*FUN(x/base) 
}

# Function to calculate Shannon Diversity Index
SHDI <- function(x, na.rm = TRUE) {
    # Convert to numeric vector if matrix
    if (is.matrix(x)) x <- as.numeric(x)
    if (na.rm) x <- na.omit(x)
    
    # Count all occurrences of each unique temperature
    props <- table(x) / length(x)
    
    # Calculate Shannon diversity index
    result <- -sum(props * log(props), na.rm = TRUE)
    
    return(result)
}

# Function to calculate Simpson Diversity Index
SIDI <- function(x, na.rm = TRUE) {
    # Convert to numeric vector if matrix
    if (is.matrix(x)) x <- as.numeric(x)
    if (na.rm) {x <- na.omit(x)}
    
    # Count all occurrences of each unique temperature
    props <- table(x) / length(x)
    
    # Calculate Simpson diversity index
    result <- 1 - sum(props * props, na.rm = TRUE)
    return(result)
}

# Calculate various uniqueness metrics separately for each patch of each species
uniqueness <- function(rgb_df, nbrs_perc){
    # Create distance matrix
    dist_mat <- 
        dist(data.frame(x = rgb_df[,"red"],
                        y = rgb_df[,"green"],
                        z = rgb_df[,"blue"])) %>% 
        as.matrix()
    row.names(dist_mat) <- rgb_df[,"species"]
    colnames(dist_mat) <- rgb_df[,"species"]
    
    # Calculate number of nearest neighbours
    N_nbrs <- floor(ncol(dist_mat) * (nbrs_perc / 100))
    
    # Sort distances for each species (row)
    dist_mat <- apply(dist_mat, MARGIN = 1, sort)
    
    # Subset to the desired number of nearest neighbours (exc self)
    dist_mat <- dist_mat[2:(N_nbrs + 1),]
    
    # Calculate species (column) means
    spp_dist <- colMeans(dist_mat)
    
    # Return as a dataframe
    spp_dist <- data.frame(species = names(spp_dist),
                           uniqueness = spp_dist,
                           stringsAsFactors = FALSE, 
                           row.names = NULL)
    return(spp_dist)
}

# Function to subset colour dataframe to desired patch and species list
sub_rgb_df <- function(spp, patch, rgb_df){
    result <- rgb_df[rgb_df[,"patch"] == patch &
                         rgb_df[,"species"] %in% spp,
                     c("species", "patch", "red", "green", "blue")]
    return(result)
}

# Calculate colour uniqueness of whole body of each species
hull_overlap <- function(dat, sp){
    # Intersect hull of focal species with that of all other species
    inter <- tryCatch(
        {
            geometry::intersectn(
                # Focal species colour data
                dat[dat[, "species"] == sp, c("red", "green", "blue")], 
                # All other species
                dat[dat[, "species"] != sp, c("red", "green", "blue")])
        },
        error = function(cond) {
            return(NA)
        })
    
    # If NA, body_uniqueness is NA and calculate pattern diversity separately
    if (all(is.na(inter))) {
        pattern_diversity <- 
            convhulln(dat[dat[, "species"] == sp, c("red", "green", "blue")], 
                      "n FA")$vol
        uniqueness_body <- NA
    }else{
        # Volume of the convex hull of the focal species is pattern diversity
        pattern_diversity <- inter$ch1$vol
        # Volume of non-overlap is a measure of whole body colour uniqueness
        uniqueness_body <- inter$ch1$vol - inter$ch$vol
    }
    
    return(list(pattern_diversity = pattern_diversity,
                uniqueness_body = uniqueness_body))
}

# Function for calculating distance to centroid from given point
dist_to_centroid <- function(centroid, pt){
    return(dist(rbind(centroid, pt), method = "euclidean"))
}

# Function for calculating distance to centroid from a list of species
dist_from_list <- 
    function(centroid, sps, traits, minsp = 4, 
             trait_vars = c("red", "green", "blue")){
        # Return NA if there are not enough species
        if (length(sps) < minsp) {
            result <- as.list(rep(NA, length(trait_vars)))
            names(result) <- trait_vars
            return(result)  
        } 
        else{
            # Subset trait matrix
            sps_traits <- traits[which(traits[, "species"] %in% sps), trait_vars]
            
            # Calculate distance to centroid for each species (row)
            result <- 
                apply(sps_traits, MARGIN = 1, 
                      FUN = dist_to_centroid, centroid = centroid)
            return(result)
        }
    }

# Function to represent 3D coordinates as one number
id3D <- function(x, y, z, multiplier){
    (2^(x * multiplier)) * (3^(y * multiplier)) * (5^(z * multiplier))
}

# Function to assign coordinates to voxels
assign_voxel <- function(coords, precision, trait_vars, multiplier = 0.01){
    # Round to given precision
    coords <- 
        apply(coords, MARGIN = 2, FUN =  function(x) 
            mround(x = x, base = precision, FUN = round))
    
    # Order
    coords <- coords[order(coords[,trait_vars[1]],
                           coords[,trait_vars[2]],
                           coords[,trait_vars[3]]),]
    
    # Assign unique cell ID
    coords <- cbind(coords, 
                    voxel = id3D(x = coords[,trait_vars[1]],
                                 y = coords[,trait_vars[2]],
                                 z = coords[,trait_vars[3]], 
                                 multiplier = multiplier))
    # Return
    return(coords)
}

# Function to count coordinates within voxels
voxel_count <- function(coords, precision, trait_vars, multiplier = 0.01){
    # Assign coordinates to voxels
    coord_voxels <-  assign_voxel(coords, precision, trait_vars, multiplier)
    
    # Count coordinates inside each voxel
    voxel_count <- table(coord_voxels[,"voxel"])
    coord_voxels <- cbind(coord_voxels[match(names(voxel_count), 
                                             coord_voxels[,"voxel"]),], 
                          voxel_count)
    # Calculate proportion of coordinates inside each voxel
    # (this allows the same metric to be used for current community, of which 
    # there is only one, and future communities, of which there are 100)
    coord_voxels <-
        as.data.frame(coord_voxels) %>%
        mutate(voxel_prop = 100 * (voxel_count / sum(voxel_count)))
    
    # Return
    return(coord_voxels)
}

# Function to average across hulls
hull_avg <- function(vertices, precision, trait_vars){
    # Simplify vertices
    simp_vertices <- lapply(vertices, voxel_count, precision = 0.5, trait_vars)
    simp_vertices <- do.call("rbind", simp_vertices)
    
    # Count number of times each voxel appears in simulations
    vertex_count <- table(simp_vertices[,"voxel"])
    
    # Restrict to vertices occurring 50 or more times
    vertex_count <- vertex_count[vertex_count >= 50]
    simp_vertices <- simp_vertices[!(duplicated(simp_vertices[,"voxel"])),]
    simp_vertices <- simp_vertices[simp_vertices[,"voxel"] %in% names(vertex_count),]
    
    # Recalculate the convex hull
    avg_hull <- convhulln_ext(simp_vertices[, trait_vars])$hull_coords
    
    # Return
    return(avg_hull)
}

SES <- function(effect, SD){
    if (is.na(effect) | is.na(SD)) {
        return(NA)
    } else if (effect == 0 & SD == 0) {
        return(0)
    }else{
        return(effect / SD)
    }
}

inf2NA <- function(x) {
    x[is.nan(x)] <- NA
    x[is.infinite(x)] <- NA
    return(x)
}

# Function to randomly select species for removal
kill_random <- function(n, colour_data){
    # Shuffle species order
    colour_data <- 
        colour_data[sample(1:nrow(colour_data), 
                           size = nrow(colour_data), 
                           replace = FALSE),]
    
    # 100 random samples of the same number of species as are traded
    random <- 
        replicate(n = 100, 
                  colour_data$species[sample(1:nrow(colour_data), 
                                             size = n, 
                                             replace = FALSE)], 
                  simplify = FALSE)
    
    # 100 random samples of the same number of species as are traded,
    # weighted by probability of extinction (higher threat = more likely to be selected)
    random_threat <- 
        replicate(n = 100, 
                  colour_data$species[
                      sample(1:nrow(colour_data), 
                             size = n, 
                             replace = FALSE, 
                             prob = 1 - (colour_data$extinction_prob))], 
                  simplify = FALSE)
    
    # Return
    return(list(random = random, random_threat = random_threat))
}

# Extension of convex hull function to return coordinates of vertices & centroid
convhulln_ext <- function(df, trait_vars){
    colmat <- as.matrix(df[, trait_vars])
    hull_results <- convhulln(colmat, "FA")
    hull <- t(hull_results$hull)
    
    # Get coordinates of hull vertices
    hull_coords <-
        data.frame(x = colmat[hull, 1],
                   y = colmat[hull, 2],
                   z = colmat[hull, 3])
    # Calculate hull centroid
    hull_centroid <- colMeans(hull_coords)
    
    # Combine results
    results <- list(hull_coords = hull_coords,
                    hull_centroid = hull_centroid,
                    hull_vol = hull_results$vol)
    # Return
    return(results)
}

# Function for calculating mean colour uniqueness from a list of species
CU_from_list <- function(sps, traits, metric){
    results <- 
        mean(traits[traits[,"species"] %in% sps, metric], 
             na.rm = TRUE)
    return(results)
}

# Function for calculating CD from a list of species
CD_from_list <- function(sps, traits, minsp = 4, trait_vars = c("red", "green", "blue")){
    # Return NA if there are not enough species
    if (length(sps) < minsp) {
        return(list(hull_coords = NA,
                    hull_centroid = NA,
                    hull_vol = NA))  
    } 
    else{
        # Subset trait matrix
        sps_traits <- traits[which(traits[, "species"] %in% sps), trait_vars]
        # Calculate colour diversity for the subset trait matrix
        return(convhulln_ext(sps_traits, trait_vars = trait_vars))
    }
}

# Internal function for CD_from_random
CD_random <- function(pool, n, traits, metric, trait_vars = c("red", "green", "blue")){
    # Random sample n species from the pool
    sp.tmp <- pool[sample(1:length(pool), n, replace = FALSE)]
    
    # Calculate CU
    CU_results <- 
        CU_from_list(sps = sp.tmp, traits = traits, metric = "CU_c_spMax")
    # Calculate CD
    CD_results <- 
        CD_from_list(sps = sp.tmp, traits = traits, trait_vars = trait_vars)$hull_vol
    
    # Return
    return(c(CU_results, CD_results))
}

# Calculates CD from a random selection of n species from the given species pool
# pool = species pool from which to sample
# n = n species to sample
# traits = trait data
# runs = times for replicate the random sample
CD_from_random <- function(pool, n, traits, runs = 100){
    if (runs < 1 ) {
        cat("ERROR: 'run' has to be higher than 1")
        stop()
    }
    # Do not calculate for less than 3 species or if species pool is empty
    if (n < 4 | length(pool) == 0) { 
        return(rep(NA, runs)) 
    }
    # Run
    results <- replicate(runs, CD_random(pool, n, traits))
    results <- matrix(results, ncol = 2, byrow = TRUE)
    colnames(results) <- c("CU", "CD")
    return(results)
}

# Shuffles extinction probabilities within communities and calculates colour 
# diversity for each random future community
CD_from_prob <- function(probs, sps, traits, runs = 1, minsp = 3, trait_vars = c("red", "green", "blue")){
    # build tree with species
    if (length(sps) < minsp) return(rep(NA,runs))
    else{
        # Subset future communities to the species that occur in this focal community
        probs <- probs[which(row.names(probs) %in% sps),]
        # Shuffle extinction probabilities
        probs <- probs[shuffle(probs),]
        
        CD.tmp <- apply(probs, MARGIN = 2, function(x)  
            
            CD_from_list(sps[which(x == 1)], traits, trait_vars = trait_vars))
        
        return(CD.tmp)
    }
}

# Function to calculate community colour diversity metrics
community_CD <- 
    function(i, 
             occ_matrix,
             colour_data,
             traded_species,
             realms,
             patch,
             trait_vars = c("red", "green", "blue")){
        # 1. Community setup ---------------------------------------------------
        
        # This community
        cell <- occ_matrix[i,]
        
        # Coordinates of focal cell (row)
        coords <- as.numeric(unlist(strsplit(row.names(cell), "_")))
        longitude <- coords[1]
        latitude <- coords[2]
        
        # Species present in focal cell
        spp_inds <- which(cell == 1)
        spp <- names(cell[spp_inds])
        spp_rich <- length(spp)
        
        # Traded species present in focal cell
        spp_trade <- spp[spp %in% traded_species]
        
        # Subset colour data
        colour_data_sub <- colour_data[colour_data$species %in% spp,]
        
        # If there are fewer than 3 species, don't run any simulations
        if (spp_rich < 4)
        {
            # Observed current community colour diversity
            CD_cc <- CD_from_list(sps = spp, traits = colour_data,
                                  trait_vars = trait_vars)$hull_vol
            CD_cc_trade <- CD_from_list(sps = spp_trade, traits = colour_data, 
                                        trait_vars = trait_vars)$hull_vol
            
            # Observed current mean community colour uniqueness
            # (mean across current species-level maximum uniqueness)
            CU_cc <- CU_from_list(sps = spp, traits = colour_data, 
                                  metric = "CU_c_spMax")
            CU_cc_trade <- CU_from_list(sps = spp_trade, traits = colour_data, 
                                        metric = "CU_c_spMax")
            
            # Observed current mean community pattern_diversity
            # (mean across current species-level convex hull volume)
            PD_cc <- CU_from_list(sps = spp, traits = colour_data, 
                                  metric = "pattern_diversity")
            PD_cc_trade <- CU_from_list(sps = spp_trade, traits = colour_data, 
                                        metric = "pattern_diversity")
            
            # Observed current mean community pattern_diversity
            # (mean across current species-level volume of non-overlap)
            BU_cc <- CU_from_list(sps = spp, traits = colour_data, 
                                  metric = "uniqueness_body")
            BU_cc_trade <- CU_from_list(sps = spp_trade, traits = colour_data, 
                                        metric = "uniqueness_body")
            
            # Bind
            results <- data.frame(
                longitude = longitude,
                latitude = latitude, 
                spN = length(spp),
                CD_cc = CD_cc,
                CD_cc_trade = CD_cc_trade,
                CD_cc_ses_glb = NA,
                CD_cc_ses_rlm = NA,
                CU_cc = CU_cc,
                CU_cc_trade = CU_cc_trade,
                CU_cc_ses_glb = NA,
                CU_cc_ses_rlm = NA,
                # _______
                # Current
                # - PATTERN DIVERSITY & BODY UNIQUENESS
                # _______
                # Community current pattern diversity
                PD_cc = PD_cc,
                PD_cc_trade = PD_cc_trade,
                # Community current body uniqueness
                BU_cc = BU_cc,
                BU_cc_trade = BU_cc_trade)
            return(results)
        } 
        else
        {
            # 2. Current colour diversity ------------------------------------------
            
            # Observed current community colour diversity
            CD_cc <- CD_from_list(sps = spp, traits = colour_data,
                                  trait_vars = trait_vars)$hull_vol
            CD_cc_trade <- CD_from_list(sps = spp_trade, traits = colour_data, 
                                        trait_vars = trait_vars)$hull_vol
            
            # Observed current mean community colour uniqueness
            # (mean across current species-level maximum uniqueness)
            CU_cc <- CU_from_list(sps = spp, traits = colour_data, 
                                  metric = "CU_c_spMax")
            CU_cc_trade <- CU_from_list(sps = spp_trade, traits = colour_data, 
                                        metric = "CU_c_spMax")
            
            # Observed current mean community pattern_diversity
            # (mean across current species-level convex hull volume)
            PD_cc <- CU_from_list(sps = spp, traits = colour_data, 
                                  metric = "pattern_diversity")
            PD_cc_trade <- CU_from_list(sps = spp_trade, traits = colour_data, 
                                        metric = "pattern_diversity")
            
            # Observed current mean community pattern_diversity
            # (mean across current species-level volume of non-overlap)
            BU_cc <- CU_from_list(sps = spp, traits = colour_data, 
                                  metric = "uniqueness_body")
            BU_cc_trade <- CU_from_list(sps = spp_trade, traits = colour_data, 
                                        metric = "uniqueness_body")
            
            # Random current community colour diversity (global species pool)
            CD_cc_glb <- 
                CD_from_random(
                    # Global species pool
                    pool = spp_colour$species, 
                    # Same species richness as current cell
                    n = spp_rich, 
                    traits = colour_data,
                    # 100 iterations
                    runs = 100)
            
            # Mean random colour diversity across iterations
            CD_cc_glb_mean <- colMeans(CD_cc_glb, na.rm = TRUE)
            # SD random colour diversity across iterations
            CD_cc_glb_sd <- apply(CD_cc_glb, MARGIN = 2, sd, na.rm = TRUE)
            # Current community SES colour diversity 
            # (observed - null mean) / null SD
            CD_cc_ses_glb <- SES(CD_cc - CD_cc_glb_mean["CD"], CD_cc_glb_sd["CD"])
            # Current community SES colour uniqueness
            CU_cc_ses_glb <- SES(CU_cc - CD_cc_glb_mean["CU"], CD_cc_glb_sd["CU"])
            
            # Random current community colour diversity (realm species pool)
            if (is.na(realms[i])) {
                CD_cc_rlm_mean <- NA
                CD_cc_rlm_sd <- NA
                CD_cc_ses_rlm <- NA
                CU_cc_ses_rlm <- NA
            } else{
                CD_cc_rlm <-
                    CD_from_random(
                        # Realm species pool
                        pool = names(which(colSums(occ_matrix[which(realms == realms[i]),], 
                                                   na.rm = TRUE) > 0)), 
                        # Same species richness as current cell
                        n = spp_rich, 
                        traits = colour_data,
                        # 100 iterations
                        runs = 100)
                CD_cc_rlm_mean <- colMeans(CD_cc_rlm, na.rm = TRUE)
                CD_cc_rlm_sd <- apply(CD_cc_rlm, MARGIN = 2, sd, na.rm = TRUE)
                CD_cc_ses_rlm <- SES(CD_cc - CD_cc_rlm_mean["CD"], CD_cc_rlm_sd["CD"])
                CU_cc_ses_rlm <- SES(CU_cc - CD_cc_rlm_mean["CU"], CD_cc_rlm_sd["CU"])
            }
            
            # 3. Gather results ----------------------------------------------------
            results <- data.frame(
                longitude = longitude,
                latitude = latitude, 
                spN = length(spp),
                # _______
                # Current
                # - COLOUR DIVERSITY
                # _______
                # Community current colour diversity
                CD_cc = CD_cc,
                CD_cc_trade = CD_cc_trade,
                # Community SES current colour diversity (global species pool)
                CD_cc_ses_glb = CD_cc_ses_glb,
                # Community SES current colour diversity (realm species pool)
                CD_cc_ses_rlm = CD_cc_ses_rlm,
                # _______
                # Current
                # - COLOUR UNIQUENESS
                # _______
                # Community current colour uniqueness
                CU_cc = CU_cc,
                CU_cc_trade = CU_cc_trade,
                # Community SES current colour uniqueness (global species pool)
                CU_cc_ses_glb = CU_cc_ses_glb,
                # Community SES current colour uniqueness (realm species pool)
                CU_cc_ses_rlm = CU_cc_ses_rlm,
                # _______
                # Current
                # - PATTERN DIVERSITY & BODY UNIQUENESS
                # _______
                # Community current pattern diversity
                PD_cc = PD_cc,
                PD_cc_trade = PD_cc_trade,
                # Community current body uniqueness
                BU_cc = BU_cc,
                BU_cc_trade = BU_cc_trade)
            
            return(results)
        }
    }   



# Function to calculate community colour diversity metrics
community_CD_future <- 
    function(i,
             patch,
             scenario,
             iter,
             kill_species,
             occ_matrix,
             colour_data,
             trait_vars = c("red", "green", "blue")){
        # 1. Community setup ---------------------------------------------------
        
        # This community
        cell <- occ_matrix[i,]
        
        # Coordinates of focal cell (row)
        coords <- as.numeric(unlist(strsplit(row.names(cell), "_")))
        longitude <- coords[1]
        latitude <- coords[2]
        
        # Species present in focal cell
        spp_inds <- which(cell == 1)
        spp <- names(cell[spp_inds])
        spp_rich <- length(spp)
        
        # Subset colour data
        colour_data_sub <- colour_data[colour_data$species %in% spp,]
        
        # If there are fewer than 3 species, don't run any simulations
        if (spp_rich < 4)
        { 
            return(NA) 
        } 
        else 
        {
            # Remove traded species
            kill_species_i <- kill_species[[scenario]][[iter]]
            kill_prop <- names(kill_species)[scenario]
            spp_traded <- spp[which(!(spp %in% kill_species_i))]
            
            # Create a random samples of the same size as the 
            # species pool after removal of traded species
            spp_random = replicate(
                n = 100,
                expr = sample(spp, 
                              size = length(spp_traded), 
                              replace = FALSE), 
                simplify = FALSE)
            
            # Observed current community colour diversity
            CD_cc <- CD_from_list(sps = spp, traits = colour_data,
                                  trait_vars = trait_vars)$hull_vol
            # Observed current mean community colour uniqueness
            # (mean across current species-level maximum uniqueness)
            CU_cc <- CU_from_list(sps = spp, traits = colour_data, 
                                  metric = "CU_c_spMax")
            
            # Future community diversity
            CD_fc <- CD_from_list(sps = spp_traded, 
                                  traits = colour_data,
                                  trait_vars = trait_vars)$hull_vol
            # Future community uniqueness (based on current uniqueness values)
            CU_fc <- CU_from_list(sps = spp_traded, 
                                  traits = colour_data, 
                                  metric = "CU_c_spMax")
            
            # Random future community diversity
            CD_frc <- 
                unlist(lapply(spp_random, function(iter)
                    CD_from_list(sps = iter, 
                                 traits = colour_data,
                                 trait_vars = trait_vars)$hull_vol))
            # Random future community uniqueness
            CU_frc <- 
                unlist(lapply(spp_random, function(iter)
                    CU_from_list(sps = iter, 
                                 traits = colour_data, 
                                 metric = "CU_c_spMax")))
            
            # Colour loss
            CD_fc_loss = CD_cc - CD_fc
            CD_frc_loss = CD_cc - CD_frc
            CU_fc_loss = CU_cc - CU_fc
            CU_frc_loss = CU_cc - CU_frc
            
            # Calculate SES
            CD_loss_ses <- 
                SES(effect = CD_fc_loss - mean(CD_frc_loss, na.rm = TRUE),
                    SD = sd(CD_frc_loss, na.rm = TRUE))
            CU_loss_ses <- 
                SES(effect = CU_fc_loss - mean(CU_frc_loss, na.rm = TRUE),
                    SD = sd(CU_frc_loss, na.rm = TRUE))
            
            return(data.frame(longitude = longitude,
                              latitude = latitude,
                              patch = patch,
                              kill_prop = kill_prop,
                              iter = iter,
                              CD_fc_loss = CD_fc_loss,
                              CD_loss_ses = CD_loss_ses,
                              CU_fc_loss = CU_fc_loss,
                              CU_loss_ses = CU_loss_ses))
            
        }
    }
