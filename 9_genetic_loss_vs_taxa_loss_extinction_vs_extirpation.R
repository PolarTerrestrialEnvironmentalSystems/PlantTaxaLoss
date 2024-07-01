# TAXALOSS Project - script to supplementary information 
# Script from Ying Liu
# Script last update - 15.04.2023

# Layout of the script: 

# 0 - test the ASV extirpated or even likely globally extinct
# you need:The species list of 7 ASVs (Supplementary table 4)

# 1 - test if the genetic loss reflect the real taxa loss or not
# you need:the outputs from previous scripts (resampled_data_final_1000iterations.rda)
# the species and corresponding ASV in Sibala_2023 database(2023_03_14_Table_before_fasta_file_final.csv)
# species_loss_artifical_timeframe (species_loss_artifical_timeframe.csv)
# the ASVs and taxa corresponding file (ASVs_taxa_after_louvain_community.csv)

###############################################################################
# 0 - test the ASV extirpated or even likely globally extinct

###############################################################################
# Remove all the objects we created so far
rm(list = ls()) 

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # state that everything is reads as character 

# Load all needed packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(rgbif)

# test for the first ASV
key <- name_suggest(q='Libanotis buchtormensis', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key, limit = 10000)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Apiaceae', rank='family')$data$key[1]

Apiaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Apiacea1<-data.frame(unique(Apiaceae_data$data$species))

###
key <- name_suggest(q=' Seseli schrenkianum', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Apiaceae', rank='family')$data$key[1]

Apiaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = 6720)

Apiacea2<-data.frame(unique(Apiaceae_data$data$species))

###
key <- name_suggest(q='Tilingia ajanensis', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Apiaceae', rank='family')$data$key[1]

Apiaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = 6720)

Apiacea3<-data.frame(unique(Apiaceae_data$data$species))

Apiacea<-rbind(Apiacea1,Apiacea2,Apiacea3)
colnames(Apiacea)<-c("Apiaceae")
###
coordinate<-data.frame(Latitude=c(55,90,90,55,55),
                       Longitude=c(150,150,50,50,150))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Apiacea <- occ_search(geometry = polygon_wkt, taxonKey = 6720)

Alaska_Apiacea1<-na.omit(data.frame(unique(Alaska_Apiacea$data$species)))


coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(180,180,150,150,180))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Apiacea <- occ_search(geometry = polygon_wkt, taxonKey = 6720)

Alaska_Apiacea2<-na.omit(data.frame(unique(Alaska_Apiacea$data$species)))

coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(-140,-140,-180,-180,-140))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Apiacea <- occ_search(geometry = polygon_wkt, taxonKey = 6720)

Alaska_Apiacea3<-na.omit(data.frame(unique(Alaska_Apiacea$data$species)))

Alaska_Apiacea<-rbind(Alaska_Apiacea1,Alaska_Apiacea2,Alaska_Apiacea3)

colnames(Alaska_Apiacea)<-c("Apiaceae")

Apiacea_samples<-data.frame(sample(Apiacea$Apiaceae,50))
colnames(Apiacea_samples)<-c("Apiaceae")

Apiacea_overlap<-data.frame(Apiacea_samples[Apiacea_samples$Apiaceae%in%Alaska_Apiacea$Apiaceae,])

colnames(Apiacea_overlap)<-c("Apiacea_overlap")

#overlap rate 20%

# test for the second ASV
key <- name_suggest(q='Arbelaezaster ellsworthii', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key, limit = 10000)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Asteraceae', rank='family')$data$key[1]

Asteraceae_data <- occ_search(geometry = polygon_wkt, taxonKey = 3065)

Asteraceae1<-na.omit(data.frame(unique(Asteraceae_data$data$species)))
colnames(Asteraceae1)<-c("Asteraceae")

###
key <- name_suggest(q='Ekmaniopappus', rank='genus')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Asteraceae', rank='family')$data$key[1]

Asteraceae_data <- occ_search(geometry = polygon_wkt, taxonKey = 3065)

Asteraceae2<-na.omit(data.frame(unique(Asteraceae_data$data$species)))

colnames(Asteraceae2)<-c("Asteraceae")


###
key <- name_suggest(q='Filago nevadensis', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Asteraceae', rank='family')$data$key[1]

Asteraceae_data <- occ_search(geometry = polygon_wkt, taxonKey = 3065)

Asteraceae3<-na.omit(data.frame(unique(Asteraceae_data$data$species)))
colnames(Asteraceae3)<-c("Asteracaeae")

###

key <- name_suggest(q='Gamochaeta lulioana', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Asteraceae', rank='family')$data$key[1]

Asteraceae_data <- occ_search(geometry = polygon_wkt, taxonKey = 3065)

Asteraceae4<-na.omit(data.frame(unique(Asteraceae_data$data$species)))
colnames(Asteraceae4)<-c("Asteracaeae")

###

key <- name_suggest(q='Pseudoclappia arenaria', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Asteraceae', rank='family')$data$key[1]

Asteraceae_data <- occ_search(geometry = polygon_wkt, taxonKey = 3065)

Asteraceae5<-na.omit(data.frame(unique(Asteraceae_data$data$species)))
colnames(Asteraceae5)<-c("Asteracaeae")
names(Asteraceae1) <- names(Asteraceae2) <- names(Asteraceae3) <- names(Asteraceae4) <- names(Asteraceae5) <- colnames(Asteraceae1)

Asteraceae<-rbind(Asteraceae1,Asteraceae2,Asteraceae3,Asteraceae4,Asteraceae5)


###
coordinate<-data.frame(Latitude=c(55,90,90,55,55),
                       Longitude=c(150,150,50,50,150))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Asteraceae <- occ_search(geometry = polygon_wkt, taxonKey = 3065)

Alaska_Asteraceae1<-na.omit(data.frame(unique(Alaska_Asteraceae$data$species)))


coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(180,180,150,150,180))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Asteraceae <- occ_search(geometry = polygon_wkt, taxonKey = 3065)

Alaska_Asteraceae2<-na.omit(data.frame(unique(Alaska_Asteraceae$data$species)))

coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(-140,-140,-180,-180,-140))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Asteraceae <- occ_search(geometry = polygon_wkt, taxonKey = 3065)

Alaska_Asteraceae3<-na.omit(data.frame(unique(Alaska_Asteraceae$data$species)))

Alaska_Asteraceae<-rbind(Alaska_Asteraceae1,Alaska_Asteraceae2,Alaska_Asteraceae3)

colnames(Alaska_Asteraceae)<-c("Asteraceaee")

Asteraceae_samples<-data.frame(sample(Asteraceae$Asteraceae,100))
colnames(Asteraceae_samples)<-c("Asteraceaee")

Alaska_Asteraceae_samples<-data.frame(sample(Alaska_Asteraceae$Asteraceaee,100))
colnames(Alaska_Asteraceae_samples)<-c("Asteraceaee")

Asteraceae_overlap<-data.frame(Asteraceae_samples[Asteraceae_samples$Asteraceaee%in%Alaska_Asteraceae_samples$Asteraceaee,])

colnames(Asteraceae_overlap)<-c("Asteraceae_overlap")
#overlap rate 2%

# test for the third and forth ASV
key <- name_suggest(q='Eritrichium caucasicum', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key, limit = 10000)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Boraginaceae', rank='family')$data$key[1]

Boraginaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Boraginaceae<-na.omit(data.frame(unique(Boraginaceae_data$data$species)))
colnames(Boraginaceae)<-c("Boraginaceae")
###

coordinate<-data.frame(Latitude=c(55,90,90,55,55),
                       Longitude=c(150,150,50,50,150))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Boraginaceae <- occ_search(geometry = polygon_wkt, taxonKey = 2498)

Alaska_Boraginaceae1<-na.omit(data.frame(unique(Alaska_Boraginaceae$data$species)))


coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(180,180,150,150,180))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Boraginaceae <- occ_search(geometry = polygon_wkt, taxonKey = 2498)

Alaska_Boraginaceae2<-na.omit(data.frame(unique(Alaska_Boraginaceae$data$species)))

coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(-140,-140,-180,-180,-140))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Boraginaceae <- occ_search(geometry = polygon_wkt, taxonKey = 2498)

Alaska_Boraginaceae3<-na.omit(data.frame(unique(Alaska_Boraginaceae$data$species)))

Alaska_Boraginaceae<-rbind(Alaska_Boraginaceae1,Alaska_Boraginaceae2,Alaska_Boraginaceae3)

colnames(Alaska_Boraginaceae)<-c("Boraginaceaee")

Boraginaceae_samples<-data.frame(sample(Boraginaceae$Boraginaceae,28))
colnames(Boraginaceae_samples)<-c("Boraginaceaee")

Alaska_Boraginaceae_samples<-data.frame(sample(Alaska_Boraginaceae$Boraginaceaee,28))
colnames(Alaska_Boraginaceae_samples)<-c("Boraginaceaee")

Boraginaceae_overlap<-data.frame(Boraginaceae_samples[Boraginaceae_samples$Boraginaceaee%in%Alaska_Boraginaceae_samples$Boraginaceaee,])

colnames(Boraginaceae_overlap)<-c("Boraginaceae_overlap")

#Boraginaceae overlap rate 25% 7/28

# test for the fifth ASV
key <- name_suggest(q='Rumex suffruticosus', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key, limit = 10000)

unique(species_distribution$data$species)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Polygonaceae', rank='family')$data$key[1]

Polygonaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Polygonaceae<-na.omit(data.frame(unique(Polygonaceae_data$data$species)))
colnames(Polygonaceae)<-c("Polygonaceae")
###

coordinate<-data.frame(Latitude=c(55,90,90,55,55),
                       Longitude=c(150,150,50,50,150))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Polygonaceae <- occ_search(geometry = polygon_wkt, taxonKey = 2416)

Alaska_Polygonaceae1<-na.omit(data.frame(unique(Alaska_Polygonaceae$data$species)))


coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(180,180,150,150,180))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Polygonaceae <- occ_search(geometry = polygon_wkt, taxonKey = 2416)

Alaska_Polygonaceae2<-na.omit(data.frame(unique(Alaska_Polygonaceae$data$species)))

coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(-140,-140,-180,-180,-140))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Polygonaceae <- occ_search(geometry = polygon_wkt, taxonKey = 2416)

Alaska_Polygonaceae3<-na.omit(data.frame(unique(Alaska_Polygonaceae$data$species)))

Alaska_Polygonaceae<-rbind(Alaska_Polygonaceae1,Alaska_Polygonaceae2,Alaska_Polygonaceae3)

colnames(Alaska_Polygonaceae)<-c("Polygonaceaee")

Polygonaceae_samples<-data.frame(sample(Polygonaceae$Polygonaceae,29))
colnames(Polygonaceae_samples)<-c("Polygonaceaee")

Alaska_Polygonaceae_samples<-data.frame(sample(Alaska_Polygonaceae$Polygonaceaee,29))
colnames(Alaska_Polygonaceae_samples)<-c("Polygonaceaee")

Polygonaceae_overlap<-data.frame(Polygonaceae_samples[Polygonaceae_samples$Polygonaceaee%in%Alaska_Polygonaceae_samples$Polygonaceaee,])

colnames(Polygonaceae_overlap)<-c("Polygonaceae_overlap")

#Polygonaceae overlap rate 27.58% 8/29


# test for the sixth ASV
key <- name_suggest(q='Adenostyles alliariae', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key, limit = 10000)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Caprifoliaceae', rank='family')$data$key[1]

Caprifoliaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = 6710)

Caprifoliaceae1<-na.omit(data.frame(unique(Caprifoliaceae_data$data$species)))
colnames(Caprifoliaceae1)<-c("Caprifoliaceae")

###
key <- name_suggest(q='Valeriana apula', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Caprifoliaceae', rank='family')$data$key[1]

Caprifoliaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Caprifoliaceae2<-na.omit(data.frame(unique(Caprifoliaceae_data$data$species)))

colnames(Caprifoliaceae2)<-c("Caprifoliaceae")


###
key <- name_suggest(q='Valeriana dioica L.', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Caprifoliaceae', rank='family')$data$key[1]

Caprifoliaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Caprifoliaceae3<-na.omit(data.frame(unique(Caprifoliaceae_data$data$species)))
colnames(Caprifoliaceae3)<-c("Caprifoliaceae")

###

key <- name_suggest(q='Valeriana pyrenaica', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Caprifoliaceae', rank='family')$data$key[1]

Caprifoliaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Caprifoliaceae4<-na.omit(data.frame(unique(Caprifoliaceae_data$data$species)))
colnames(Caprifoliaceae4)<-c("Caprifoliaceae")

###

key <- name_suggest(q='Valeriana saliunca', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

key <- name_suggest(q='Caprifoliaceae', rank='family')$data$key[1]

Caprifoliaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Caprifoliaceae5<-na.omit(data.frame(unique(Caprifoliaceae_data$data$species)))
colnames(Caprifoliaceae5)<-c("Caprifoliaceae")
names(Caprifoliaceae1) <- names(Caprifoliaceae2) <- names(Caprifoliaceae3) <- names(Caprifoliaceae4) <- names(Caprifoliaceae5) <- colnames(Caprifoliaceae1)

Caprifoliaceae<-rbind(Caprifoliaceae1,Caprifoliaceae2,Caprifoliaceae3,Caprifoliaceae4,Caprifoliaceae5)


###
coordinate<-data.frame(Latitude=c(55,90,90,55,55),
                       Longitude=c(150,150,50,50,150))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Caprifoliaceae <- occ_search(geometry = polygon_wkt, taxonKey = 6710)

Alaska_Caprifoliaceae1<-na.omit(data.frame(unique(Alaska_Caprifoliaceae$data$species)))


coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(180,180,150,150,180))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Caprifoliaceae <- occ_search(geometry = polygon_wkt, taxonKey = 6710)

Alaska_Caprifoliaceae2<-na.omit(data.frame(unique(Alaska_Caprifoliaceae$data$species)))

coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(-140,-140,-180,-180,-140))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Caprifoliaceae <- occ_search(geometry = polygon_wkt, taxonKey = 6710)

Alaska_Caprifoliaceae3<-na.omit(data.frame(unique(Alaska_Caprifoliaceae$data$species)))

Alaska_Caprifoliaceae<-rbind(Alaska_Caprifoliaceae1,Alaska_Caprifoliaceae2,Alaska_Caprifoliaceae3)

colnames(Alaska_Caprifoliaceae)<-c("Caprifoliaceaee")

Caprifoliaceae_samples<-data.frame(sample(Caprifoliaceae$Caprifoliaceae,29))
colnames(Caprifoliaceae_samples)<-c("Caprifoliaceaee")

Alaska_Caprifoliaceae_samples<-data.frame(sample(Alaska_Caprifoliaceae$Caprifoliaceaee,29))
colnames(Alaska_Caprifoliaceae_samples)<-c("Caprifoliaceaee")

Caprifoliaceae_overlap<-data.frame(Caprifoliaceae_samples[Caprifoliaceae_samples$Caprifoliaceaee%in%Alaska_Caprifoliaceae_samples$Caprifoliaceaee,])

colnames(Caprifoliaceae_overlap)<-c("Caprifoliaceae_overlap")
#Caprifoliaceae overlap rate 24%  7/29



# test for the seventh ASV
key <- name_suggest(q='Potentilla ancistrifolia', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key, limit = 10000)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")


key <- name_suggest(q='Rosaceae', rank='family')$data$key[1]

Rosaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Rosaceae1<-na.omit(data.frame(unique(Rosaceae_data$data$species)))
colnames(Rosaceae1)<-c("Rosaceae")
###

key <- name_suggest(q='Potentilla dickinsii', rank='species')$data$key[1]

species_distribution <- occ_search(taxonKey = key, limit = 10000)

coordinate<-na.omit(data.frame(Latitude=species_distribution$data$decimalLatitude,
                               Longitude=species_distribution$data$decimalLongitude))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")


key <- name_suggest(q='Rosaceae', rank='family')$data$key[1]

Rosaceae_data <- occ_search(geometry = polygon_wkt, taxonKey = key)

Rosaceae2<-na.omit(data.frame(unique(Rosaceae_data$data$species)))
colnames(Rosaceae2)<-c("Rosaceae")

Rosaceae<-rbind(Rosaceae1,Rosaceae2)

###
coordinate<-data.frame(Latitude=c(55,90,90,55,55),
                       Longitude=c(150,150,50,50,150))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Rosaceae <- occ_search(geometry = polygon_wkt, taxonKey = 5015)

Alaska_Rosaceae1<-na.omit(data.frame(unique(Alaska_Rosaceae$data$species)))


coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(180,180,150,150,180))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Rosaceae <- occ_search(geometry = polygon_wkt, taxonKey = 5015)

Alaska_Rosaceae2<-na.omit(data.frame(unique(Alaska_Rosaceae$data$species)))

coordinate<-data.frame(Latitude=c(40,90,90,40,40),
                       Longitude=c(-140,-140,-180,-180,-140))

BHL = chull(coordinate[, c("Latitude", "Longitude")]) 
BHL = c(BHL, BHL[1])

convex_hull_points <- coordinate[BHL, c("Longitude", "Latitude")]
convex_hull_points <- rbind(convex_hull_points, convex_hull_points[1, ])

polygon_wkt <- paste("POLYGON((", paste(convex_hull_points$Longitude, 
                                        convex_hull_points$Latitude, 
                                        collapse = ","), "))", sep = "")

Alaska_Rosaceae <- occ_search(geometry = polygon_wkt, taxonKey = 5015)

Alaska_Rosaceae3<-na.omit(data.frame(unique(Alaska_Rosaceae$data$species)))

Alaska_Rosaceae<-rbind(Alaska_Rosaceae1,Alaska_Rosaceae2,Alaska_Rosaceae3)

colnames(Alaska_Rosaceae)<-c("Rosaceaee")

Rosaceae_samples<-data.frame(sample(Rosaceae$Rosaceae,100))
colnames(Rosaceae_samples)<-c("Rosaceaee")

Alaska_Rosaceae_samples<-data.frame(sample(Alaska_Rosaceae$Rosaceaee,100))
colnames(Alaska_Rosaceae_samples)<-c("Rosaceaee")

Rosaceae_overlap<-data.frame(Rosaceae_samples[Rosaceae_samples$Rosaceaee%in%Alaska_Rosaceae_samples$Rosaceaee,])

colnames(Rosaceae_overlap)<-c("Rosaceae_overlap")

#Rosaceae overlap rate 4% 4/100


###############################################################################
# 1 - test if the genetic loss reflect the real taxa loss or not
###############################################################################
rm(list = ls())

# Import the resampled data
load("resampled_data_final_1000iterations.rda")
# Import the database information
database<-read.csv("2023_03_14_Table_before_fasta_file_final.csv")

ASVs_taxa<-read.csv("ASVs_taxa_after_louvain_community.csv")%>%
  select(uniq_label,NUC_SEQ)

resampleAll <- smpl_raw
look <- resampleAll[[1]][[1]]

# Set age timeslice: 
age <- tibble(timeslice = c(13:0), age = c(1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000, 17000, 19000, 21000, 23000, 25000, 27000))
df_long_for_rar <- read_delim("df_long_for_rar.csv", delim = ",", col_names = T) %>%
  subset(read_counts > 0)

df_long_for_rar %>% #subset(mean_age < 28001) %>% 
  dplyr::select(mean_age, sample_id) %>%
  #summarise(sum = sum(read_counts)) %>%
  distinct() 

taxa_info <- read_delim("final_dataset_before_resampling_family_merge_new.csv", delim = ",", col_names = T) %>%
  dplyr::select(identity, family, genus, species, scientific_name, group, nb_in_group) %>% mutate(uniq = 1:nrow(.)) %>% 
  mutate(level = ifelse(is.na(species), ifelse(is.na(genus), "family", "genus"), "species")) %>%
  mutate(genus = ifelse(is.na(genus), family, genus)) %>% mutate(species = ifelse(is.na(species), paste(genus, uniq, sep = " "), paste(species, uniq, sep = " "))) %>% mutate(identity = ifelse(identity == 0.9, "cand", "1")) %>%
  mutate(key = paste(nb_in_group, identity, family, scientific_name, sep = "_"))

###############################################################################
# 1 - Make the abundance table
###############################################################################
abund_table <- lapply(resampleAll, function(x) x[[1]])
abund_table <- lapply(abund_table, function(x) replace(x, is.na(x), 0))

# Change the info on NcountTab adn save as different outputs
NcountTab <- lapply(abund_table, function(x) x %>% pivot_longer(`13`:`0`) %>% subset(value > 0) %>% 
                      mutate(name = as.numeric(name)) %>% setNames(c("uniq_label", "time_slice", "N")) %>% 
                      arrange(time_slice) %>% as.data.table() %>%
                      dcast(uniq_label ~ time_slice, value.var = "N")  %>%
                      mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
                      subset(grepl("_cand_", uniq_label)) %>% # to set if we want to calculate the extinction rate
                      #subset(grepl("_1_", uniq_label)) %>% # to set if we want to calculate the extirpation rate
                      right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label"))
a<-NcountTab[[1]]

NcountTab1 <- lapply(NcountTab, function(x) {
  x$uniq_label <- gsub("grass_|herb_|tree_|shrub_", "", x$uniq_label)
  merge(x, ASVs_taxa, by="uniq_label")
})

b<-NcountTab1[[1]]


NcountTab2 <- lapply(NcountTab1, function(x) {
  x$uniq_label <- x$NUC_SEQ
  x%>%select(-NUC_SEQ)
})

c<-NcountTab2[[1]]

nSAMPLE <- 1

reapTab = list()
for (i in 1:nSAMPLE) {
  reapTab[[i]] <- rbind(
    apply(t(apply(NcountTab2[[i]][,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab2[[i]][,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab2[[i]][,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab2[[i]][,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
}

# got the ASV lost in every timeslice
save(reapTab,file="reapTab_extinction.Rdata")

# using the 13 artifical timeslcie (interval between 14 timeslices), and use 
# the ASV and species from the Sibala database, to select the ASVs randomly to 
# set the ASV lost (the number of lost ASVs in every artifical timeslice is similar
# with the ASV lost from the corresponding timeslice), then for every artifical 
# timeslice, the number of lost species (corresponding to the lost ASV in the
# artifical timeslice) were checked; then in every artifical timeslice,the observed
# species lost and expected species lost were compared, and Proportion of loss 
# species relative to the expectation base on simulation were got.

species_all<-list()

nresample<-1000

for (r in 1:nresample) {
  print(r)
  
  species<-read.csv("species.csv",check.names = FALSE)
  rownames(species)<-species$Age
  species<-species[,-c(1)]
  
  sequence<-as.data.frame(unique(database$DNA_SEQ))
  colnames(sequence)<-"sequence"
  
  removed_sequences <- vector("list", 13)
  
  reapTabi<-reapTab[[r]]
  new_row<-t(data.frame(
    row1 = c("T1", "T2", "T3", "T4", "T5", "T6"," T7", "T8", "T9", "T10", "T11", "T12", "T13"),
    row2= c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1,0)
  ))
  
  reapTabi<-rbind(new_row,reapTabi)
  rownames(reapTabi)<-c("Timeslice",
                        "potential number of timeslice to recover",
                        "total disappear",
                        "reappear never",
                        "reappear after 13 timeslices",
                        "reappear after 12 timeslices",
                        "reappear after 11 timeslices",
                        "reappear after 10 timeslices",
                        "reappear after 9 timeslices",
                        "reappear after 8 timeslices",
                        "reappear after 7 timeslices",
                        "reappear after 6 timeslices",
                        "reappear after 5 timeslices",
                        "reappear after 4 timeslices",
                        "reappear after 3 timeslices",
                        "reappear after 2 timeslices",
                        "reappear after 1 timeslices")
  reapTabi<-as.data.frame(reapTabi)
  colnames(reapTabi)<-colnames(test)[c(2:14)]
  num_slices <- 13
    for (i in 1:num_slices) {
    removed_sequences[[i]] <- as.data.frame(sequence[sample(1:nrow(sequence), reapTabi[3, i]),]) 
    
    colnames(removed_sequences[[i]]) <- "sequence"
    
    sequence <- as.data.frame(sequence[!sequence$sequence %in% removed_sequences[[i]]$sequence,])
    
    colnames(sequence) <- "sequence"
    
    remove_taxa <- as.data.frame(database[database$DNA_SEQ %in% removed_sequences[[i]]$sequence,])
    
    non_remove_taxa<-as.data.frame(database[database$DNA_SEQ %in% sequence$sequence,])
    
    remove_taxa <- as.data.frame(remove_taxa[!remove_taxa$ncbi_species %in% non_remove_taxa$ncbi_species,])
    
    species[3, i] <- as.numeric(length(unique(remove_taxa$ncbi_species)))
    
    if (i > 1) {
      for (j in 1:(i-1)) {
        back_data <- as.data.frame(removed_sequences[[j]][sample(1:nrow(removed_sequences[[j]]), reapTabi[(18-i+j), j]),])
        
        colnames(back_data) <- "sequence"
        
        back_taxa <- as.data.frame(database[database$DNA_SEQ %in% back_data$sequence,])
        
        non_remove_taxa1<-as.data.frame(database[database$DNA_SEQ %in% sequence$sequence,])
        
        back_taxa<- as.data.frame(back_taxa[!back_taxa$ncbi_species %in% non_remove_taxa1$ncbi_species,])
        
        species[18-i+j, j] <- as.numeric(length(unique(back_taxa$ncbi_species)))
        
        sequence <- rbind(sequence, back_data)
        removed_sequences[[j]] <- as.data.frame(removed_sequences[[j]][!removed_sequences[[j]]$sequence %in% back_data$sequence,])
        colnames(removed_sequences[[j]]) <- "sequence"
      }
    }
  }
  
  
  species[-4,][species[-4,] == ""] <- 0
  
  species<-species[-c(1),]
  
  species<-type.convert(species,as.is=TRUE)
  
  for (col in 1:ncol(species)) {
    
    species[3, col] <- species[2, col] - sum(species[4:16, col])
  }
  
  species[species<0]<-0
  
  species_observe_reappear<- as.data.frame(100*(species[c(2),]-species[c(3),])/(species[c(2),]))
  
  species_expected_reappear<- t(data.frame(100*sum(species[16,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[15,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[14,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[13,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[12,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[11,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[10,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[9,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[8,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[7,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[6,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[5,c(1:4)])/sum(species[2,c(1:4)]),
                                           100*sum(species[4,c(1:4)])/sum(species[2,c(1:4)])))
  
  cumulative_species <- species_expected_reappear[, 1]
  
  for (i in 2:nrow(species_expected_reappear)) {
    cumulative_species[i] <- cumulative_species[i-1] + species_expected_reappear[i, 1]
  }
  
  cumulative_species<-as.data.frame(cumulative_species)
  
  colnames(cumulative_species)<-"expected reappearance"
  
  cumulative_species$observe<-t(species_observe_reappear[,c(13:1)])  
  
  cumulative_species<-as.data.frame(cumulative_species)
  
  cumulative_species$difference<-cumulative_species$`expected reappearance`-cumulative_species$observe
  cumulative_species$round<-r
  cumulative_species$age<-rep(seq(1000,25000,by=2000))
  species_all[[r]]<-cumulative_species
}

species_all1 <- lapply(species_all, function(df) {
  colnames(df) <- c("expected","observed","difference")
  return(df)
})

df_species_9timeslice<-do.call(rbind, species_all)
colnames(df_species_9timeslice) <- c("expected","observed","difference","round","age")

write.csv(df_species_9timeslice,"extinction_df_species_9timeslice.csv")


summary_data<-df_species_9timeslice%>%
  filter(age!=1000&age!=25000)%>%
  filter(expected<100)%>%
  setNames(c("expected", "observed","difference","round","age")) %>%
  group_by(age) %>%
  summarise(median = median(difference,na.rm = TRUE),
            q95 = quantile(difference, 0.95,na.rm = TRUE),
            q5 = quantile(difference, 0.05,na.rm = TRUE))

#make the plot
plot1<-ggplot(summary_data) +
  geom_line(aes(x = age, y = median/100)) +
  geom_ribbon(aes(x = age,ymin = q95/100, ymax = q5/100), alpha = 0.2) +
  labs(x = "Age / year BP", y = "Possibility of species loss") +
  geom_point(aes(x = age,y=median/100), size=3) +
  geom_hline(yintercept=0, linetype = "dashed") +
  scale_x_reverse(breaks=seq(0, 25000, 2000)) +
  xlab("Age (cal. yrs BP)") +
  ylab("Proportion of loss species relative to the expectation base on simulation") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot1

ggsave(plot1,filename="extinction_9timeslice.pdf",width = 12,height=8)


































