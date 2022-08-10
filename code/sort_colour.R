
# Setup -------------------------------------------------------------------

library(dplyr)

# Define function to get hex code 
hex <- function(r, g, b) rgb(r / 255, g / 255, b / 255)

# Read in raw colour data
colour <- read.csv("data/raw/colour_raw.csv")

# Add HSV & colour categories ---------------------------------------------

colour <-
    colour %>% 
    mutate(
        # Identify which sex is more colourful
        colourful_sex = 
            case_when(Male_plumage_score - Female_plumage_score >= 0 ~ "male",
                      Male_plumage_score - Female_plumage_score < 0 ~ "female"),
        # Match to hexadecimal code
        hex = hex(red, green, blue))

# Calculate HSV
colour[, c("hue", "saturation", "value")] <- 
    t(rgb2hsv(r = colour$red, g = colour$green, b = colour$blue))
colour <-
    colour %>% 
    mutate(
        # Rescale
        hue = hue * 360,
        saturation = saturation * 100,
        value = value * 100,
        # Assign to colour categories
        colour_cat = 
            case_when(value <= 20 ~ "dark",
                      value > 20 & saturation <= 10 ~ "light",
                      value > 20 & saturation > 10 &
                          hue >= 345 | (hue >= 0 & hue < 15) ~ "red",
                      # Special category for brown vs. orange
                      # (orange only if both saturation AND value > 90)
                      hue >= 15 & hue < 45 &
                          value >= 90 & saturation >= 90 ~ "orange",
                      # (brown if either saturation OR value < 90)
                      hue >= 15 & hue < 45 &
                          value > 20 & saturation > 10 &
                          (value < 90 | saturation < 90) ~ "brown",
                      value > 20 & saturation > 10 &
                          hue >= 45 & hue < 75 ~ "yellow",
                      value > 20 & saturation > 10 &
                          hue >= 75 & hue < 105 ~ "chartreuse_green",
                      value > 20 & saturation > 10 &
                          hue >= 105 & hue < 135 ~ "green",
                      value > 20 & saturation > 10 &
                          hue >= 135 & hue < 165 ~ "spring_green",
                      value > 20 & saturation > 10 &
                          hue >= 165 & hue < 195 ~ "cyan",
                      value > 20 & saturation > 10 &
                          hue >= 195 & hue < 225 ~ "azure",
                      value > 20 & saturation > 10 &
                          hue >= 225 & hue < 255 ~ "blue",
                      value > 20 & saturation > 10 &
                          hue >= 255 & hue < 285 ~ "violet",
                      value > 20 & saturation > 10 &
                          hue >= 285 & hue < 315 ~ "magenta",
                      value > 20 & saturation > 10 &
                          hue >= 315 & hue < 345 ~ "rose"))

# Relevel
colour <- 
    colour %>% 
    mutate(patch = factor(patch,
                          levels = c("nape",
                                     "crown",
                                     "forehead",
                                     "throat",
                                     "upper_breast",
                                     "lower_breast")),
           colour_cat = factor(colour_cat,
                               levels = c("light",
                                          "dark",
                                          "brown",
                                          "red",
                                          "orange",
                                          "yellow",
                                          "chartreuse_green",
                                          "green",
                                          "spring_green",
                                          "cyan",
                                          "azure",
                                          "blue",
                                          "violet",
                                          "magenta",
                                          "rose")),
           # Broader IUCN classification
           iucn_broad = 
               case_when(iucn %in% c("LC", "NT") ~ "low_risk",
                         iucn %in% c("VU", "EN", "CR", 
                                     # Possibly extinct in the wild
                                     "CR (PEW)", 
                                     # Possibly extinct
                                     "CR (PE)") ~ "high_risk"),
           iucn = factor(iucn,
                         levels = c("LC", "NT", 
                                    "VU", "EN", "CR", 
                                    # Possibly extinct in the wild
                                    "CR (PEW)", 
                                    # Possibly extinct
                                    "CR (PE)", 
                                    "EW", "EX", "DD"),
                         labels = c("LC", "NT", 
                                    "VU", "EN", 
                                    rep("CR",3), 
                                    "EW", "EX", "DD")),
           traded = factor(traded,
                           levels = c(0, 1),
                           labels = c("no", "yes")))

# Write -------------------------------------------------------------------
saveRDS(colour, file = "data/raw/colour_sorted.Rds")
