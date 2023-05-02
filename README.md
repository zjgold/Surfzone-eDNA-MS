# Surfzone-eDNA-MS
A Comparison of Biomonitoring Methodologies for Surf Zone Fish Communities

Zachary Gold1 ¶*, McKenzie Q. Koch1¶, Nicholas K. Schooler2, Kyle A. Emery2, Jenifer E. Dugan2, Robert J. Miller2, Henry M. Page2, Donna M. Schroeder3, David M. Hubbard2, Jessica R. Madden2, Stephen G. Whitaker2,4, and Paul H. Barber1

1Department of Ecology and Evolutionary Biology, University of California Los Angeles, CA, USA

2Marine Science Institute, University of California Santa Barbara, CA, USA 

3Bureau of Ocean Energy Management, Camarillo, CA, USA

4Channel Islands National Park, Ventura, CA, USA


*  Corresponding Author
Email: zjgold@ucla.edu 
Phone: +1(310)-795-0020 

¶ These authors contributed equally to this work.

## Abstract
Surf zones are highly dynamic marine ecosystems that are subject to increasing anthropogenic and climatic pressures, posing multiple challenges for biomonitoring. Traditional methods such as seines and hook and line surveys are often labor intensive, taxonomically biased, and can be physically hazardous. Emerging techniques, such as baited remote underwater video (BRUV) and environmental DNA (eDNA) are promising nondestructive tools for assessing marine biodiversity in surf zones of sandy beaches. Here we compare the relative performance of beach seines, BRUV, and eDNA in characterizing community composition of bony (teleost) and cartilaginous (elasmobranch) fishes of surf zones at 18 open coast sandy beaches in southern California. Seine and BRUV surveys captured overlapping, but distinct fish communities with 50% (18/36) of detected species shared. BRUV surveys more frequently detected larger species (e.g. sharks and rays) while seines more frequently detected one of the most abundant species, barred surfperch (Amphistichus argenteus). In contrast, eDNA metabarcoding captured 88.9% (32/36) of all fishes observed in seine and BRUV surveys plus 57 additional species, including 15 that frequent surf zone habitats. On average, eDNA detected over 5 times more species than BRUVs and 8 times more species than seine surveys at a given site. eDNA approaches also showed significantly higher sensitivity than seine and BRUV methods and more consistently detected 31 of the 32 (96.9%) jointly observed species across beaches. The four species detected by BRUV/seines, but not eDNA were only resolved at higher taxonomic ranks (e.g. Embiotocidae surfperches and Sygnathidae pipefishes). In frequent co-detection of species between methods limited comparisons of richness and abundance estimates, highlighting the challenge of comparing biomonitoring approaches. Despite potential for improvement, results overall demonstrate that eDNA can provide a cost-effective tool for long-term surf zone monitoring that complements data from seine and BRUV surveys, allowing more comprehensive surveys of vertebrate diversity in surf zone habitats.

**Preprint:** [*Link*](https://www.biorxiv.org/content/10.1101/2021.11.19.469341v1.full)
**Accepted at PLOS ONE**

### Description of Important Files
1. 20221231_ucsb_2018_ms_analysis.Rmd - Code to conduct analysis
2. 20230103_SOM.R - Code to conduct site occupancy modeling
3. figure_1_code.Rmd - Code to generate Figure 1 map
4. decontam/ - Directory of decontamination processing of Anacapa Output tables.
4.a 20230102_generalized_decontam_ucsb_2018_exploration_tech_reps.R - Code to run decontamination
5. data/ - Directory of data
5.a. anacapa_20230101/ - Directory of Anacapa Toolkit output ASV tables

All other files are intermediate data tables used across analyses. Appologies in advance for the mess. At least everything you would need to recreate the analyses from scratch are here.
