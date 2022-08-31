# install packages


## note at lines 395, LDhap API token = "************" 
## must be requested from https://ldlink.nci.nih.gov/?tab=apiaccess
## and inputted in place of ************


install.packages("LDlinkR") # import package to query the LDhap API

install.packages("dplyr") 

install.packages("assertr")
install.packages("readr")
install.packages("httr")
install.packages('comprehenr')



#call libraries

library(LDlinkR)
library(purrr) # to flatten list
library(assertr)
library(readr)
library(httr)
library(comprehenr) # for list comprehension; to_list



## CREATE A LIST OF MAJOR POPULATION GROUPS
cross_continent <- list(all = "ALL", # all populations
                        european = "EUR", # european
                        african = "AFR", # african 
                        admixed_american = "AMR", # ad mixed american 
                        east_asian = "EAS", # east asian 
                        south_asian = "SAS" # south asian 
                       )


## CREATE A LIST OF SUB-POPULATION GROUPS IN EACH MAJOR POPULATION GROUP

## EUROPEAN 
european <- list(utah_NW_europe = "CEU", # utah residents from north and west europe
                toscani = "TSI", # toscani om italia
                finnish = "FIN", # finnish in england
                british = "GBR", # british in england and scotland
                iberian = "IBS" # iberian population in spain
                )

## INTRA-AFRICA
ad_mixed_american <- list(mexican = "MXL", # mexican ancestry from los angeles usa
                         puerto_ricans = "PUR", # Puerto Ricans from Puerto Rico
                         colombians = "CLM", # colombians from medellin colombia
                         peruvians = "PEL" # peruvians from lima peru
                         )

## EAST ASIAN
east_asian <- list(han_chinese = "CHB", # han chinese in beijing china
                  japanese = "JPT", # japanese in tokyo japan
                  southern_han = "CHS", # southern han chinese
                  chinese_dai = "CDX", # chinese dai in xishuangbanna china
                  kinh = "KHV" # kinh in ho chi minh city vietnam
                  )

## SOUTH ASIAN
south_asian <- list(gujarati = "GIH", # gujarati indian from houston texas
                   punjabi = "PJL", # punjabi from lahore pakistan
                   bengali = "BEB", # bengali from bangladesh
                   sri_lankan_tamil = "STU", # sri lankan tamil from the uk
                   indian_telugu = "ITU" # indian telugu from the uk
                   )



## pharmGKB data has been put in a github repo
## pull the data from github

alprazolam_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/alprazolam_Clin.tsv"))
alprazolam_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/alprazolam_Var.tsv"))
amitriptyline_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/amitriptyline_Clin.tsv"))
amitriptyline_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/amitriptyline_Var.tsv"))
carbamazepine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/carbamazepine_Clin.tsv"))
carbamazepine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/carbamazepine_Var.tsv"))
chlorpromazine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/chlorpromazine_Clin.tsv"))
chlorpromazine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/chlorpromazine_Var.tsv"))
citalopram_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/citalopram_Clin.tsv"))
citalopram_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/citalopram_Var.tsv"))
clozapine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/clozapine_Clin.tsv"))
clozapine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/clozapine_Var.tsv"))
diazepam_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/diazepam_Clin.tsv"))
diazepam_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/diazepam_Var.tsv"))
fluoxetine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/fluoxetine_Clin.tsv"))
fluoxetine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/fluoxetine_Var.tsv"))
fluphenazine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/fluphenazine_Clin.tsv"))
fluphenazine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/fluphenazine_Var.tsv"))
haloperidol_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/haloperidol_Clin.tsv"))
haloperidol_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/haloperidol_Var.tsv"))
imipramine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/imipramine_Clin.tsv"))
imipramine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/imipramine_Var.tsv"))
lithium_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/lithium_Clin.tsv"))
lithium_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/lithium_Var.tsv"))
lorazepam_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/lorazepam_Clin.tsv"))
lorazepam_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/lorazepam_Var.tsv"))
midazolam_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/midazolam_Clin.tsv"))
midazolam_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/midazolam_Var.tsv"))
olanzapine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/olanzapine_Clin.tsv"))
olanzapine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/olanzapine_Var.tsv"))
paroxetine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/paroxetine_Clin.tsv"))
paroxetine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/paroxetine_Var.tsv"))
phenobarbital_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/phenobarbital_Clin.tsv"))
phenobarbital_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/phenobarbital_Var.tsv"))
phenytoin_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/phenytoin_Clin.tsv"))
phenytoin_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/phenytoin_Var.tsv"))
quetiapine_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/quetiapine_Clin.tsv"))
quetiapine_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/quetiapine_Var.tsv"))
risperidone_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/risperidone_Clin.tsv"))
risperidone_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/risperidone_Var.tsv"))
sertraline_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/sertraline_Clin.tsv"))
sertraline_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/sertraline_Var.tsv"))
valproic_acid_Clin <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/valproic_acid_Clin.tsv"))
valproic_acid_Var <- read_tsv(url("https://raw.github.com/BBLinus/pharmgenomics/main/Clin_var_annotation/valproic_acid_Var.tsv"))



## list of dataframe containing the drugs

pharm_tables <- list(alprazolam_Clin = alprazolam_Clin,
                     alprazolam_Var = alprazolam_Var, 
                     amitriptyline_Clin = amitriptyline_Clin,
                     amitriptyline_Var = amitriptyline_Var,
                     carbamazepine_Clin = carbamazepine_Clin,
                     carbamazepine_Var = carbamazepine_Var, 
                     chlorpromazine_Clin = chlorpromazine_Clin, 
                     chlorpromazine_Var = chlorpromazine_Var, 
                     citalopram_Clin = citalopram_Clin, 
                     citalopram_Var = citalopram_Var, 
                     clozapine_Clin = clozapine_Clin, 
                     clozapine_Var = clozapine_Var, 
                     diazepam_Clin = diazepam_Clin, 
                     diazepam_Var = diazepam_Var, 
                     fluoxetine_Clin = fluoxetine_Clin, 
                     fluoxetine_Var = fluoxetine_Var, 
                     fluphenazine_Clin = fluphenazine_Clin, 
                     fluphenazine_Var = fluphenazine_Var, 
                     haloperidol_Clin = haloperidol_Clin, 
                     haloperidol_Var = haloperidol_Var, 
                     imipramine_Clin = imipramine_Clin, 
                     imipramine_Var = imipramine_Var, 
                     lithium_Clin = lithium_Clin,
                     lithium_Var = lithium_Var,
                     lorazepam_Clin = lorazepam_Clin,
                     lorazepam_Var = lorazepam_Var,
                     midazolam_Clin = midazolam_Clin,
                     midazolam_Var = midazolam_Var,
                     olanzapine_Clin = olanzapine_Clin,
                     olanzapine_Var = olanzapine_Var,
                     paroxetine_Clin = paroxetine_Clin,
                     paroxetine_Var = paroxetine_Var,
                     phenobarbital_Clin = phenobarbital_Clin,
                     phenobarbital_Var = phenobarbital_Var,
                     phenytoin_Clin = phenytoin_Clin,
                     phenytoin_Var = phenytoin_Var,
                     quetiapine_Clin = quetiapine_Clin,
                     quetiapine_Var = quetiapine_Var,
                     risperidone_Clin = risperidone_Clin,
                     risperidone_Var = risperidone_Var,
                     sertraline_Clin = sertraline_Clin,
                     sertraline_Var = sertraline_Var,
                     valproic_acid_Clin = valproic_acid_Clin,
                     valproic_acid_Var = valproic_acid_Var)

length(pharm_tables)
names(pharm_tables)


# CLINICAL ANNOTATIONS, VARIANT ANNOTATIONS EXCLUDED
clin_annot <- to_list(for (i in seq(1, length(pharm_tables), by=2)) {
    print(pharm_tables[[i]])
    }
                      )

# clinical annotation names
clin_annot_names <- to_list(for (i in seq(1, length(names(pharm_tables)), by=2)) {
    print(names(pharm_tables)[[i]])
    }
                            )

names(clin_annot) <- c(clin_annot_names)
names(clin_annot)

# number of drugs
length(clin_annot)

## CLINICAL ANNOTATIONS WITH ONLY MODERATE TO HIGH LEVEL OF EVIDENCE (1A, 1B, 2A, 2B)

alprazolam <- clin_annot$alprazolam_Clin[clin_annot$alprazolam_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
amitriptyline <- clin_annot$amitriptyline_Clin[clin_annot$amitriptyline_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
carbamazepine <- clin_annot$carbamazepine_Clin[clin_annot$carbamazepine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
chlorpromazine <- clin_annot$chlorpromazine_Clin[clin_annot$chlorpromazine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
citalopram <- clin_annot$citalopram_Clin[clin_annot$citalopram_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
clozapine <- clin_annot$clozapine_Clin[clin_annot$clozapine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
diazepam <- clin_annot$diazepam_Clin[clin_annot$diazepam_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
fluoxetine <- clin_annot$fluoxetine_Clin[clin_annot$fluoxetine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
fluphenazine <- clin_annot$fluphenazine_Clin[clin_annot$fluphenazine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
haloperidol <- clin_annot$haloperidol_Clin[clin_annot$haloperidol_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
imipramine <- clin_annot$imipramine_Clin[clin_annot$imipramine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
lithium <- clin_annot$lithium_Clin[clin_annot$lithium_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
lorazepam <- clin_annot$lorazepam_Clin[clin_annot$lorazepam_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
midazolam <- clin_annot$midazolam_Clin[clin_annot$midazolam_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
olanzapine <- clin_annot$olanzapine_Clin[clin_annot$olanzapine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
paroxetine <- clin_annot$paroxetine_Clin[clin_annot$paroxetine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
phenobarbital <- clin_annot$phenobarbital_Clin[clin_annot$phenobarbital_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
phenytoin <- clin_annot$phenytoin_Clin[clin_annot$phenytoin_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
quetiapine <- clin_annot$quetiapine_Clin[clin_annot$quetiapine_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
risperidone <- clin_annot$risperidone_Clin[clin_annot$risperidone_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
sertraline <- clin_annot$sertraline_Clin[clin_annot$sertraline_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]
valproic_acid <- clin_annot$valproic_acid_Clin[clin_annot$valproic_acid_Clin$Level %in% c('1A', '1B', '2A', '2B'), ]


clinical_evidence <- list(alprazolam = alprazolam,
                          amitriptyline = amitriptyline,
                          carbamazepine = carbamazepine,
                          chlorpromazine = chlorpromazine,
                          citalopram = citalopram,
                          clozapine = clozapine,
                          diazepam = diazepam, 
                          fluoxetine = fluoxetine, 
                          fluphenazine = fluphenazine,
                          haloperidol = haloperidol,
                          imipramine = imipramine,
                          lithium = lithium, 
                          lorazepam = lorazepam,
                          midazolam = midazolam, 
                          olanzapine = olanzapine,
                          paroxetine = paroxetine,
                          phenobarbital = phenobarbital,
                          phenytoin = phenytoin, 
                          quetiapine = quetiapine,
                          risperidone = risperidone,
                          sertraline = sertraline,
                          valproic_acid = valproic_acid) 

clinical_evidence

alprazolam
amitriptyline
carbamazepine
chlorpromazine
citalopram
clozapine
diazepam
fluoxetine
fluphenazine
haloperidol
imipramine
lithium
lorazepam
midazolam
olanzapine
paroxetine
phenobarbital
phenytoin
quetiapine
risperidone
sertraline
valproic_acid



## LIST OF IDENTIFIED VARIANTS

variants <- c(alprazolam$Variant, amitriptyline$Variant, carbamazepine$Variant, chlorpromazine$Variant, citalopram$Variant,
              clozapine$Variant, diazepam$Variant, fluoxetine$Variant, fluphenazine$Variant, haloperidol$Variant, imipramine$Variant,
              lithium$Variant, lorazepam$Variant, midazolam$Variant, olanzapine$Variant, paroxetine$Variant, phenobarbital$Variant,
              phenytoin$Variant, quetiapine$Variant, risperidone$Variant, sertraline$Variant, valproic_acid$Variant)

length(list_variants)

variants_list <- unique(list_variants)
variants_list

# number of variants
length(variants_list)


## Tag SNPs for Each Haplotype
tag_list <- list('CYP2D6*1' = '', # normal function gene
                 'CYP2D6*1xN' = '', # copy number variation
                 'CYP2D6*2' = c('rs1058164', 'rs16947', 'rs1135840'), # C>S(G/C), G>A, C>G
                 'CYP2D6*2xN' = '', # copy number variation
                 'CYP2D6*3' = c('rs35742686'), # T>delT # ignore rs1135824 (T>Y(C/T) from pharmGKB keep rs35742686
                 'CYP2D6*4' = c('rs1065852', 'rs28371703', 'rs28371704', 'rs1058164', 'rs3892097', 'rs16947', 'rs1058172',
                                 'rs1135833', 'rs1135835', 'rs28371735', 'rs766507177', 'rs74478221', 'rs75467367', 'rs1135840'), 
                 # G>R(G/A), G>K(G/T), T>Y(C/T), C>S(G/C), C>T, G>R(G/A), C>Y(C/T), G>S(G/C), T>Y(C/T), G>R(G/A), T>K(G/T), C>Y(C/T), G>S(G/C), C>S(G/C)
                 'CYP2D6*5' = '', # Whole gene deletion
                 'CYP2D6*6' = c('rs5030655', 'rs1135840'), # A>delA, C>S(G/C)
                 'CYP2D6*9' = 'rs5030656', # TCT>delTCT                 
                 'CYP2D6*10' = c('rs1065852', 'rs1058164', 'rs1135840'), # G>A, C>G, C>G
                 'CYP2D6*17' = c('rs28371706', 'rs1058164', 'rs16947', 'rs1135840'), # G>A, C>G, G>A, C>G
                 'CYP2D6*14' = c('rs1058164', 'rs5030865', 'rs16947', 'rs1135840'), # C>G, C>T, G>A, C>G
                 'CYP2D6*41' = c('rs1058164', 'rs16947', 'rs28371725', 'rs1135840'), # C>G, G>A, C>T, C>G
                 
                 'CYP2C19*1' = 'rs3758581', # A>G
                 'CYP2C19*2' = c('rs12769205', 'rs58973490', 'rs4244285', 'rs3758581'), # A>G, G>R(G/A), G>A, A>G
                 'CYP2C19*3' = c('rs4986893', 'rs3758581'), # G>A, A>G
                 'CYP2C19*4' = c('rs12248560', 'rs28399504', 'rs3758581'), # C>Y(C/T), A>G, A>G
                 'CYP2C19*17' = c('rs12248560', 'rs3758581'), # C>T, A>G
                 
                 'CYP2C9*1' = '', # normal function gene
                 'CYP2C9*2' = 'rs1799853', # C>T
                 'CYP2C9*3' = 'rs1057910', # A>C
                 'CYP2C9*5' = 'rs28371686', # C>G
                 'CYP2C9*6' = 'rs9332131', # A>delA
                 'CYP2C9*8' = 'rs7900194', # G>A
                 'CYP2C9*11' = 'rs28371685', # C>T
                 'CYP2C9*13' = 'rs72558187', # T>C
                 'CYP2C9*14' = 'rs72558189', # G>A
                 'CYP2C9*16' = 'rs72558192', # A>G
                 'CYP2C9*29' = 'rs182132442', # C>A
                 'CYP2C9*31' = 'rs57505750', # T>C
                 'CYP2C9*33' = 'rs200183364', # G>A
                 'CYP2C9*37' = 'rs564813580', # A>G
                 'CYP2C9*39' = 'rs762239445', # G>T
                 'CYP2C9*42' = 'rs12414460', # G>A
                 'CYP2C9*43' = 'rs767576260', # C>T 
                 'CYP2C9*45' = 'rs199523631', # C>T
                 'CYP2C9*50' = '', # variants not in dbSNP
                 'CYP2C9*52' = 'rs988617574', # C>G
                 'CYP2C9*55' = 'rs1250577724', # C>A
                 
                 'CYP3A4*1' = '', # normal function gene 
                 'CYP3A4*20' = 'rs67666821', #T>TT 
                 'CYP3A4*22' = 'rs35599367', #G>A
                 
                 # HLA-B*15:02:01
                 # Cross-ethnicity tagging SNPs for HLA alleles associated with adverse drug reaction
                 'HLA-B*15:02_cross' = 'rs10484555', # C 
                 # A high-resolution HLA and SNP haplotype map for disease association studies in the extended human MHC
                 'HLA-B*15:02_chb' = c('rs3909184', 'rs2844682'), # G>Y(C/T), G>A # 2,4 # no yri, no ceu, no jpt 
                 
                 # Therefore, the HLA-B*15:02 allele is defined by its sequence rather than single coding or protein variations. 
                 # If there is strong linkage disequilibrium between one or more SNPs and a specific HLA allele, 
                 # the presence of these SNPs (tag SNPs) may be used for HLA typing (51).
                 
                 'rs3812718-T' = 'rs3812718', # C>T
                 
                 # HLA-B*40:01:01
                 # A high-resolution HLA and SNP haplotype map for disease association studies in the extended human MHC
                 'HLA-B*40:01_ceu' = c('rs4711240', 'rs1265110'), # T>M(A/C), C>T # 3,1 # no yri
                 'HLA-B*40:01_chb' = c('rs2844580', 'rs2523694'), # T>C, G>A, # 3,2
                 'HLA-B*40:01_jpt' = 'rs7381331', # ALSO rs2523567 # C>G # 3
                 
                 # Predicting HLA alleles from high-resolution SNP data in three Southeast Asian populations
                 'HLA-B*40:01_chs' = 'rs2523605', # (T) CHS: Southern Han Chinese
                 'HLA-B*40:01_mas_ins' = c('rs2523605', 'rs2844569', 'rs2507984'), # (T, T, A) MAS: Southeast Asian Malays, INS: South Asian Tamil Indians
                 
                 # HLA-A*31:01:02
                 # A high-resolution HLA and SNP haplotype map for disease association studies in the extended human MHC
                 'HLA-A*31:01_ceu' = 'rs1061235', # A>T # 4 # europe ceu, no yoruba yri
                 'HLA-A*31:01_chb' = c('rs1061235', 'rs3823318'), # A>T, C>G # 2,1 chinese chb
                 'HLA-A*31:01_jpt' = c('rs1061235', 'rs1150739'), # A>T, T>M(A/C) # 1,1 japanese jpt
                 
                 # Predicting HLA alleles from high-resolution SNP data in three Southeast Asian populations
                 'HLA-A*31:01_chs' = c('rs2523979', 'rs1150740'),  #  (A, A) CHS: Southern Han Chinese
                 'HLA-A*31:01_ins' = 'rs3909134', # (G) INS: South Asian Tamil Indians  
                 
                 # Genome-wide association study identifies HLA-A*3101 allele as a genetic risk factor for 
                 # carbamazepine-induced cutaneous adverse drug reactions in Japanese population 
                 # https://academic.oup.com/hmg/article/20/5/1034/2385803
                 'HLA-A*31:01_jap' = 'rs16333021', # A/G # Japanese
                 
                 # HLA-B*15:11:01
                 # A high-resolution HLA and SNP haplotype map for disease association studies in the extended human MHC
                 'HLA-B*15:11_chb' = 'rs2074491' # T>C # 4 # no europe ceu, no yoruba yri, no japan jpt
                 )

length(tag_list)

# RECALL POPULATION LIST
cross_continent
european
intra_africa
ad_mixed_american
east_asian
south_asian


## Identify each SNP

 
snpDF <- function (snps, pop) {
    df <- LDhap(snps = snps,
                     pop = pop, 
                     token = "************") # using the LDhap API token 
    # which can be requested at https://ldlink.nci.nih.gov/?tab=apiaccess
    return(cbind(pop,df))
}


## IDENTIFY THE VARIANTS IN THE 5 MAJOR POPULATION GROUPS
for (j in tag_list) {
    for (i in cross_continent) {
        try({
            df <- snpDF(snps=j, pop=i)
            print (df)
            write.table(df, file='cross_continent_SNP.txt', 
                        append=TRUE, sep=',', col.names=TRUE, row.names=FALSE)         
        }, silent=TRUE)
    }
}



## VARIANTS IN THE EUROPEAN SUBPOPULATIONS
for (j in tag_list) {
    for (i in european) {
        try({
            df <- snpDF(snps=j, pop=i)
            print (df)
            write.table(df, file='european_SNP.txt', 
                        append=TRUE, sep=',', col.names=TRUE, row.names=FALSE)         
        }, silent=TRUE)
    }
}



## VARIANTS IN THE AFRICAN SUBPOPULATIONS
for (j in tag_list) {
    for (i in intra_africa) {
        try({
            df <- snpDF(snps=j, pop=i)
            print (df)
            write.table(df, file='intra_africa_SNP.txt', 
                        append=TRUE, sep=',', col.names=TRUE, row.names=FALSE)         
        }, silent=TRUE)
    }
}



## VARIANTS IN THE ADMIXED AMERICAN SUBPOPULATIONS
for (j in tag_list) {
    for (i in ad_mixed_american) {
        try({
            df <- snpDF(snps=j, pop=i)
            print (df)
            write.table(df, file='ad_mixed_american_SNP.txt', 
                        append=TRUE, sep=',', col.names=TRUE, row.names=FALSE)         
        }, silent=TRUE)
    }
}


## VARIANTS IN THE EAST ASIAN SUBPOPULATIONS
for (j in tag_list) {
    for (i in east_asian) {
        try({
            df <- snpDF(snps=j, pop=i)
            print (df)
            write.table(df, file='east_asian_SNP.txt', 
                        append=TRUE, sep=',', col.names=TRUE, row.names=FALSE)         
        }, silent=TRUE)
    }
}



## VARIANTS IN THE SOUTH ASIAN SUBPOPULATIONS
for (j in tag_list) {
    for (i in south_asian) {
        try({
            df <- snpDF(snps=j, pop=i)
            print (df)
            write.table(df, file='south_asian_SNP.txt', 
                        append=TRUE, sep=',', col.names=TRUE, row.names=FALSE)         
        }, silent=TRUE)
    }
}


print('all done')
