# Abundance_HarvestedPops

Code and data for: 
Keever, A. C., J. D. Kelly, and B. S. Cohen. Estimating abundance of harvested populations at the management unit scale. _In preparation_

This repository provides all code and raw data (and associated metadata) to recreate the analyses for the above manuscript (Keever et al., _In preparation_). 

This readme file was last updated on 2023-06-05 by Allison C Keever

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
  Current Institution: 
	Address: 

	Co-Author
	Name: Bradley S. Cohen
	Institution: Tennessee Technological University
	Address: College of Arts and Sciences, Tennessee Tech University, 1 William L Jones Dr., Cookeville TN 38505, USA

3. Date of data collection: 2005-2019 

4. Geographic location of data collection: Tennessee, United States

5. Information about funding sources that supported the collection of the data: 
	This project was supported by Wildlife Restoration Grants administered by the U.S. Fish and Wildlife Service, Wildlife and Sport Fish Restoration Program: Partnering to fund 
  conservation and connect people with nature. 

## Metadata: Data & File Overview

1. Description of dataset
	These data were collected by Tennessee Wildlife Resources Agency from hunter check stations, online reporting by hunters, and from mail-in surveys sent to license holders in Tennessee.  

2. File List: 
	TN_AbundanceEstimation.R: R script to estimate abundance of deer in Tennessee
  	SimulatedDeer_Analyses.R: R script to simulate deer populations, monitoring data, and estimates of abundance
  	TN_IPM_Base.txt: Model text file for estimating abundance of deer in Tennessee 
  	IPM_HarvestRecon_Base.txt: Model text file for estimating abundance of simulated populations
  	Rscripts/Functions_AbundanceEstimation.R: R script that includes functions to estimate abundance of deer in Tennessee 
  	Rscripts/SimPopFunctions.R: R script that includes functions for simulating deer populations
  	Rscripts/TN_harvest_data_prep.R: R script that preps raw harvest data from Tennessee for analysis
  	Data/AF_Harvest_Summary.csv: Summary of age-at-harvest data for female deer
  	Data/Am_Harvest_Summary.csv: Summay of age-at-harvest data for male deer
  	Data/dkmaster.RData: Raw deer harvest data from hunter check stations and reported harvest in Tennessee
  	Data/dmus.csv: Deer management unit counties to relate county of harvest to deer management units
  	Data/HarvestEstimates.csv: Estimated harvest from social surveys sent out to license holders
  	Data/Reported_Harvest_Summary.csv: Summary of reported harvest
  	Results/ParameterEstimates.csv: A complete list of estimates of abundance and demographic rates for deer in Tennessee
  
## Metadata: Data-specific information for - AF_Harvest_Summary.csv

1. Number of variables: 6

2. Number of cases/rows: 135

3. Variable List: 
	DMU: the deer management unit the harvest occured in; character
	Year: the harvest year (not calendar year) the harvest occured in; integer - range = 2005-2019
	Fawn: the number of fawns harvested; integer
	Yearling: the number of yearlings harvested; integer
	Adult: the number of adults harvested; integer
	Total: the total number of females harvested; integer

4. Missing data codes: None

5. Specialized formats or other abbreviations used: None
  
## Metadata: Data-specific information for - Am_Harvest_Summary.csv

1. Number of variables: 6

2. Number of cases/rows: 135

3. Variable List: 
	DMU: the deer management unit the harvest occured in; character
	Year: the harvest year (not calendar year) the harvest occured in; integer - range = 2005-2019
	Fawn: the number of fawns harvested; integer
	Yearling: the number of yearlings harvested; integer
	Adult: the number of adults harvested; integer
	Total: the total number of males harvested; integer

4. Missing data codes: None

5. Specialized formats or other abbreviations used: None
  
## Metadata: Data-specific information for - dmus.csv

1. Number of variables: 6

2. Number of cases/rows: 95

3. Variable List: 
	FID: relic number code that has no meaning; integer - only value is -1
	NAME: the name of the county in Tennessee; character
	label: the full name of the designated deer management unit; character
	alias: the abbreviated name of the designated deer management unit; character
	division: the grand division in Tennessee; character
	pool; the abbreviated name to pool deer management units with similar characteristics; character

4. Missing data codes: None

5. Specialized formats or other abbreviations used: None
  
## Metadata: Data-specific information for - HarvestEstimates.csv

1. Number of variables: 6

2. Number of cases/rows: 300

3. Variable List: 
	Year: the harvest year (not calendar year) for the estimated harvest; integer - range = 2005-2019
	DMU: the deer management unit for the estimated harvest; character
	Bag: bag type (antlerless or antlered) for the estimated harvest; character
	Harv_est; estimated harvest; integer
	SD_harv_est: the standard deviation of the estimated harvest; numeric
	CL: the 95% confidence limit of the estimated harvest; integer

4. Missing data codes: NA

5. Specialized formats or other abbreviations used: State under DMU represents the estimate for the entire state of Tennessee, which includes those that had unknown locations of harvest
  
## Metadata: Data-specific information for - Reported_Harvest_Summary.csv

1. Number of variables: 4

2. Number of cases/rows: 135

3. Variable List: 
	DMU: the deer management unit the harvest occured in; character
	Year: the harvest year (not calendar year) the harvest occured in; integer - range = 2005-2019
	Antlerless: the number of antleress deer reported as harvested; integer
	Antlered: the number of antlered deer reported as harvested; integer

4. Missing data codes: None

5. Specialized formats or other abbreviations used: None
  
## Metadata: Data-specific information for - ParameterEstimates.csv

1. Number of variables: 15

2. Number of cases/rows: 4,581

3. Variable List: 
	.variable: the demographic rate the estimate is for. Hest = estimated harvest; HS = harvest survival rate; NS = natural survival rate; Report = reporting rate; N = abundance; S0 = neonate survival rate; lambda.dmu = population growth rate; rec.rate.dmu = recruitment rate
	DMU: the deer management unit the harvest occured in; character
	Year: calendar year for estimates; integer - range = 2005-2019
	Age: the age class for the estimates. F = fawn, Y = yearling; A = adult; character
	Sex: the sex class for the estimates. F = female; M = male; character
	a: number denoting age class. 1 = fawn; 2 = yearling; 3 = adult; integer
	k: number denoting year of harvest from 1 = 2005 to 15 = 2019; integer
	s: binary number denoting sex class. 1 = female; 2 = male; integer
	i: number denoting the deer management unit; integer
	.value: median estimate; numeric
	.lower: lower 90% credible interval for the estimate; numeric
	.upper: upper 90% credible interval for the estimate; numeric
	.width: width of the credible interval for the estimate (i.e., 0.9); numeric
	.point: point estimate type (i.e., median); character
	.interval: interval type (i.e., qi); character

4. Missing data codes: NA

5. Specialized formats or other abbreviations used: None
  
## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
