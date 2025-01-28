# Differentiation-and-homogenization-of-European-communities
Data and code associated with the paper (under revision) "Recovering European river invertebrate communities homogenize or differentiate depending on anthropogenic stress"
1. Code associated with analyses and figures
2. Data1_communityBasins: contains community data at the basin scale, environmental variables, and basin descriptors
   basin: basin code
   country: country name
   basin_name: basin name
   distance_mean: mean distance among sites in a basin
   distance_min: minimum distance among sites in a basin
   taxon_resolution: level of taxonomic resolution provided for all sites in a basin
   sample_units: units of abundance for the sampling sites in a basin
   number_sites: number of sites in a basin
   year: sampling year
   ts_length: length of the time series (last year - first year)
   taxa: total number of taxa across sites in a basin in a year
   total_taxa: total number of taxa across sites and years in a basin
   urban_mean: mean proportion of urban areas in the sites' upstream areas
   crop_mean: mean proportion of croplands in the sites' upstream areas
   forest_mean: mean proportion of forested areas in the sites' upstream areas
   eqr: mean ecological quality ratio (EQR) across sites in a basin in a year
   temperature_trend: mean temperature slope across sites in a basin over the time span
   lat_center: latitude of the midpoint of the basin
   lon_center: longitude of the midpoint of the basin
   taxonomic_bray: beta diversity based on Bray Curtis index and abundance-weighted taxonomic composition
   biological_bray: beta diversity based on Bray Curtis index and abundance-weighted composition of biological traits
   ecological_bray: beta diversity based on Bray Curtis index and abundance-weighted composition of ecological traits
   taxonomic_chiSq: beta diversity based on BRay Curts index and abundance-weighted taxonomic composition (Chi-Square transformed abundances)
3. Data2_communitySites: contains community data at the site scale
   basin: basin code (as in Data1)
   site: site unique ID
   country: country name
   basin_name: basin name
   latitude: site latitude
   longitude: site longitude
   year: sampling year
   abundance: total abundance
   taxon_richness: total number of distinctive taxa
   biological_richness: volume of multidimensional space occupied by biological traits
   ecological_richness: volume of multidimensional space occupied by ecological traits
   taxonomic_hill1: Shannon diversity index (or Hill number 1) represents the number of taxa weighted by the geometric mean of their proportional abundances
   biological_hill1: Shannon diversity index (or Hill number 1) represents the number of biological traits weighted by the geometric mean of their proportional abundances
   ecological_hill1: Shannon diversity index (or Hill number 1) represents the number of ecological traits weighted by the geometric mean of their proportional abundances
   eqr_site: ecological quality ratio (EQR) for a site and a year
4. Data3_ecopart: contains community data at the basin scale
   basin: basin code (as in Data1)
   country: country name
   year: sampling year
   delta_taxonomic_allYears: year-to-year change in taxonomic beta diversity
   delta_biological_allYears: year-to-year change in biological traits beta diversity
   delta_ecological_allYears: year-to-year change in ecological traits beta diversity
   taxa_addition_allYears: year-to-year change in taxonomic beta diversity due to additions
   taxa_subtraction_allYears: year-to-year change in taxonomic beta diversity due to subtractions
   biological_addition_allYears: year-to-year change in biological traits beta diversity due to additions
   biological_subtraction_allYears: year-to-year change in biological traits beta diversity due to subtractions
   ecological_addition_allYears: year-to-year change in ecological traits beta diversity due to additions
   ecological_subtraction_allYears: year-to-year change in ecological traits beta diversity due to subtractions
   delta_taxonomic_toFirst: change in taxonomic beta diversity between every year and the first year
   delta_biological_toFirst: change in biological traits beta diversity between every year and the first year
   delta_ecological_toFirst: change in ecological traits beta diversity between every year and the first year
   taxa_addition_toFirst: change in taxonomic beta diversity due to additions between every year and the first year
   taxa_subtraction_toFirst: change in taxonomic beta diversity due to subtractions between every year and the first year
   biological_addition_toFirst: change in biological traits beta diversity due to additions between every year and the first year
   biological_subtraction_toFirst: change in biological traits beta diversity due to subtractions between every year and the first year
   ecological_addition_toFirst: change in ecological traits beta diversity due to additions between every year and the first year
   ecological_subtraction_toFirst: change in ecological traits beta diversity due to subtractions between every year and the first year
   delta_taxonomic_fromLast: last-to-first year change in taxonomic beta diversity
   delta_biological_fromLast: last-to-first year change in biological traits beta diversity
   delta_ecological_fromLast: last-to-first year change in ecological traits beta diversity
   taxa_addition_fromLast: last-to-first year change in taxonomic beta diversity due to additions
   taxa_subtraction_fromLast: last-to-first year change in taxonomic beta diversity due to subtractions
   biological_addition_fromLast: last-to-first year change in biological traits beta diversity due to additions
   biological_subtraction_fromLast: last-to-first year change in biological traits beta diversity due to subtractions
   ecological_addition_fromLast: last-to-first year change in ecological traits beta diversity due to additions
   ecological_subtraction_fromLast: last-to-first year change in ecological traits beta diversity due to subtractions
5. Data4_ecopart_biotraits: contains trait data at the basin scale
   basin: basin code (as in Data1)
   year: sampling year
   biotrait: biological trait modality
   additive: change in biological trait beta diversity due to the trait addition
   subtractive: change in biological trait beta diversity due to the trait subtraction
6. Data5_ecopart_ecotraits: contains trait data at the basin scale
   basin: basin code (as in Data1)
   year: sampling year
   ecotrait: ecological trait modality
   additive: change in ecological trait beta diversity due to the trait addition
   subtractive: change in ecological trait beta diversity due to the trait subtraction
