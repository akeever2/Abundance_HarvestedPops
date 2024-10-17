[![DOI](https://zenodo.org/badge/645453351.svg)](https://zenodo.org/badge/latestdoi/645453351)

# Abundance_HarvestedPops

Code and data for: 
Keever, A. C., J. D. Kelly, G. B. Clevinger, and B. S. Cohen. Estimating abundance of harvested populations at the management unit scale. _In preparation_

This repository provides all code and raw data (and associated metadata) to recreate the analyses for the above manuscript (Keever et al., _In preparation_). Below is the metadata for files in this repository. 

This readme file was last updated on 2024-10-19 by Allison C Keever

## Metadata: General Information

1. Title of Repository: Abundance_HarvestedPops - Code and data for: Estimating abundance of harvested populations at the management unit scale.

2. Author Information  
	Author/Principal Investigator Information  
	Name: Allison C Keever  
	ORCID: 0000-0002-5194-3987  
	Institution: Tennessee Technological University  
	Address: College of Arts and Sciences, Tennessee Tech University, 1 William L Jones Dr., Cookeville TN 38505, USA  
	Email: akeever@tntech.edu  

	Co-Author  
	Name: James D. Kelly  
	Institution at time of project: Tennessee Wildlife Resources Agency  
	Address: Ellington Agricultural Center, 5105 Edmonson Pike, Nashville, TN 37211, USA  
	Current Institution:  Florida Fish and Wildlife Conservation Commission, Gainesville, Florida, United States of America

	Co-Author  
	Name: Garrett B. Clevinger  
	Institution: Tennessee Wildlife Resources Agency  
	Address: Ellington Agricultural Center, 5105 Edmonson Pike, Nashville, TN 37211, USA   

	Co-Author  
	Name: Bradley S. Cohen  
	Institution: Tennessee Technological University  
	Address: College of Arts and Sciences, Tennessee Tech University, 1 William L Jones Dr., Cookeville TN 38505, USA  

3. Date of data collection: 2005-2023 

4. Geographic location of data collection: Tennessee, United States

5. Information about funding sources that supported the collection of the data:  
	This project was supported by Wildlife Restoration Grants administered by the U.S. Fish and Wildlife Service, Wildlife and Sport Fish Restoration Program: Partnering to fund 
  conservation and connect people with nature. 

## Metadata: Data & File Overview

1. Description of dataset  
	These data were collected by Tennessee Wildlife Resources Agency from hunter check stations, online reporting by hunters, and from mail-in surveys sent to license holders in Tennessee.  

2. File List:  
	TN_AbundanceEstimation.R: R script to estimate abundance of deer in Tennessee    
  	TN_IPM_RecRate.txt: Model text file for estimating abundance of deer in Tennessee   
  	Rscripts/Functions_AbundanceEstimation.R: R script that includes functions to estimate abundance of deer in Tennessee 
  	Rscripts/TN_harvest_data_prep.R: R script that preps raw harvest data from Tennessee for analysis  
  	Data/AF_Harvest_Summary.csv: Summary of age-at-harvest data for female deer  
  	Data/Am_Harvest_Summary.csv: Summary of age-at-harvest data for male deer  
  	Data/dkmaster.RData: Raw deer harvest data from hunter check stations and reported harvest in Tennessee  
  	Data/dmus.csv: Deer management unit counties to relate county of harvest to deer management units  
  	Data/HarvestEstimates.csv: Estimated harvest from social surveys sent out to license holders  
  	Data/Reported_Harvest_Summary.csv: Summary of reported harvest    
  
## Metadata: Data-specific information for - AF_Harvest_Summary.csv

1. Number of variables: 6

2. Number of cases/rows: 114

3. Variable List:  
	DMU: the deer management unit the harvest occurred in; character  
	Year: the harvest year (not calendar year) the harvest occurred in; integer - range = 2005-2023  
	Fawn: the number of fawns harvested; integer  
	Yearling: the number of yearlings harvested; integer  
	Adult: the number of adults harvested; integer  
	Total: the total number of females harvested; integer  

4. Missing data codes: None

5. Specialized formats or other abbreviations used: None
  
## Metadata: Data-specific information for - Am_Harvest_Summary.csv

1. Number of variables: 6

2. Number of cases/rows: 114

3. Variable List:  
	DMU: the deer management unit the harvest occurred in; character  
	Year: the harvest year (not calendar year) the harvest occurred in; integer - range = 2005-2023  
	Fawn: the number of fawns harvested; integer  
	Yearling: the number of yearlings harvested; integer  
	Adult: the number of adults harvested; integer  
	Total: the total number of males harvested; integer  

4. Missing data codes: None

5. Specialized formats or other abbreviations used: None
  
## Metadata: Data-specific information for - dmus.csv

1. Number of variables: 2

2. Number of cases/rows: 95

3. Variable List:  
	DMU: The deer management unit the county resides in; character  
	County: the name of the county in Tennessee; character   

4. Missing data codes: None

5. Specialized formats or other abbreviations used: None
  
## Metadata: Data-specific information for - HarvestEstimates.csv

1. Number of variables: 6

2. Number of cases/rows: 266

3. Variable List:  
	year: the harvest year (not calendar year) for the estimated harvest; integer - range = 2005-2023  
	DMU: the deer management unit for the estimated harvest; character  
	Bag: bag type (antlerless or antlered) for the estimated harvest; character  
	Harv_est; estimated harvest; integer  
	SD_harv_est: the standard deviation of the estimated harvest; numeric  
	CL: the 95% confidence limit of the estimated harvest; integer  

4. Missing data codes: NA

5. Specialized formats or other abbreviations used: State under DMU represents the estimate for the entire state of Tennessee, which includes those that had unknown locations of harvest
  
## Metadata: Data-specific information for - Reported_Harvest_Summary.csv

1. Number of variables: 4

2. Number of cases/rows: 121

3. Variable List:  
	DMU: the deer management unit the harvest occurred in; character  
	Year: the harvest year (not calendar year) the harvest occurred in; integer - range = 2005-2019  
	Antlerless: the number of antlerless deer reported as harvested; integer  
	Antlered: the number of antlered deer reported as harvested; integer  

4. Missing data codes: None

5. Specialized formats or other abbreviations used: None
  
## Metadata: Data-specific information for - dkmaster.RData

1. Number of variables: 12

2. Number of cases/rows: 3,041,144

3. Variable List:  
	county: the Tennessee county in which the harvest occurred; character
	age: the age of the deer (0.5, 1.5, 2.5, 3.5, 4.5, and 5.5+); character
	sex: sex of the harvested deer; character
	huntyear: year of harvest based on the hunting season year, not calendar year; factor
	bag: type of bag (antlered or antlerless); character 
	HarvestDate: date of harvest; date
	AntlerPoints: antler points of the harvested deer; numeric
	Weight: weight of the harvested deer; numeric
	InsideSpreadLenghtInches: inside spread of the harvested deer in inches; numeric
	beamlength: beam length of the harvested deer; numeric
	beamcirc: beam circumference of the harvested deer; numeric
	spread: antler spread of the harvested deer; numeric

5. Missing data codes: NA

6. Specialized formats or other abbreviations used: None
  
## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
