# Differentiation-and-homogenization-of-European-communities
Data and code associated with the paper (under revision) "Recovering European river invertebrate communities homogenize or differentiate depending on anthropogenic stress"
1. Code associated with analyses and figures
2. Data 1 - Beta diversity by basin and year: contains beta-diversity data for each basin and sampling year
   basin: basin code
   year: sampling year
   tax_bray: beta diversity based on Bray Curtis index and abundance-weighted taxonomic composition
   bio_bray: beta diversity based on Bray Curtis index and abundance-weighted composition of biological traits
   eco_bray: beta diversity based on Bray Curtis index and abundance-weighted composition of ecological traits
   
3. Data 2 - Diversity indices by site and year: contains community data for each site and sampling year
   basin: basin code (as in Data1)
   site: site unique ID
   year: sampling year
   tax.richness: total number of distinctive taxa
   bio.richness: volume of multidimensional space occupied by biological traits
   eco.richness: volume of multidimensional space occupied by ecological traits
   
4. Data 3 - Environmental data per sites and year: contains environmental variables reflecting anthropogenic stress for each site and sampling year
   basin: basin code (as in Data1)
   year: sampling year
   site: site unique ID (as in Data2)
   country: country name
   latitude: decimal degrees geographic coordinate WGS84 
   longitude: decimal degrees geographic coordinate WGS84 
   eqr: ecological quality ratio
   urb.full: proportion of urban areas in the upstream area
   crop.full: proportion of agricultural areas in the upstream area
   tree.full: proportion of forested areas in the upstream area
   temp.change: slope estimate from temperature trends over the sampling period
   
5. Data 4 - Ecopart: contains additive and subtractive dynamic components of temporal beta diversity change
   basin: basin code (as in Data1)
   year: sampling year
   tax_subtraction: subtractive component of taxonomic beta diversity change
   tax_addition: additive component of taxonomic beta diversity change
   bio_subtraction: subtractive component of biological trait beta diversity change
   bio_addition: additive component of biological trait beta diversity change
   eco_subtraction: subtractive component of ecological trait beta diversity change
   eco_addition: additive component of ecological trait beta diversity change
   
7. Data 5 - Ecopart Biologica traits: contains additive and subtractive dynamic components of temporal beta diversity change calculated for each biological trait
   basin: basin code (as in Data1)
   year: sampling year
   biotrait: biological trait modality
   additive: change in biological trait beta diversity due to trait addition
   subtractive: change in biological trait beta diversity due to trait subtraction
   
9. Data 6 - Ecopart Ecologica traits: contains additive and subtractive dynamic components of temporal beta diversity change calculated for each biological trait
   basin: basin code (as in Data1)
   year: sampling year
   ecotrait: ecological trait modality
   additive: change in ecological trait beta diversity due to trait addition
   subtractive: change in ecological trait beta diversity due to trait subtraction
