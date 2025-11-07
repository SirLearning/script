main.nf: entrance of workflow
- genotype: workflows for genotype analysis
	- Assess: assess the vcf file and check the content of variant and filter others
	- Kinship: count kinship  between samples and combine different samples with similiar kinship to one taxa
	- PS: compute population structure in different methods
	- Stats: perform dimention reduce and other stats methods to get statistic values of genotypes for other analysis. Use python or r codes in this project.
- Phenotype: workflows for phenotype analysis
	- Process: get phenotype data from database and other functions, can also upload information to database, with connection through java.
	- Stats: perform statistics analysis to phenotype data to get some statistic values like BLUP and do works like combine different phenotypes for multi-traits GWAS analysis.
- Run: run GWAS analysis in different GWAS model and get result plots.