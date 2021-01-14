These are the instruction to reproduce the results of the paper called 
“Integrative analyses of TEDDY Omics data reveal lipid metabolism abnormalities, increased intracellular ROS and heightened inflammation prior to autoimmunity for type 1 diabetes”

Leandro Balzano-Nogueira, Ricardo Ramirez, Tatyana Zamkovaya, Jordan Dailey, Alexandria N. Ardissone, Srikar Chamala, Michael J. Haller, Åke Lernmark, William A. Hagopian, Marian J. Rewers, Jin-Xiong She, Jorma Toppari, Anette-G. Ziegler, Beena Akolkar, Jeffrey P. Krischer, Patrick Concannon, Mark A. Atkinson, Desmond A. Schatz, Eric W. Triplett, Ana Conesa1, and the TEDDY Study Group.


To be able to obtain the same results there are some important aspect to take into account:

1.- You will need a High performance computing center since the analyses are highly memory demanding.

2.- You must run the scripts as the following order to obtain the arrays and tables to continue obtaining the results presented in the paper:

	2.1.- VIP_NPLSDAGeneExp.R					VIP_NPLSDA analysis of the Gene Expression Data
	2.2.- VIP_NPLSDAMetabolomics.R					VIP_NPLSDA analysis of the Metabolomics Data
	2.3.- VIP_NPLSDADietBiom.R					VIP_NPLSDA analysis of the Dietary Biomarkers Data
	2.4.- VIP_NPLSDAGEMETABDB.R					VIP_NPLSDA analysis of all selected fatures to determine their descriptive capability
	2.5.- Permutationtest.R						Test to evaluate if the model retains information related to the disease
	2.6.- EnrichmentAnalysis.R					Enrichment Analysis of the individuals used, performed evaluating Cases vs Controls
	2.7.- EnrichmentAnalysisGADAIAA.R				Enrichment Analysis of the individuals used, performed evaluating GADA-first Cases vs IAA-first Cases
	2.8.- PartialCorrelationNetworks.R				Partial correlation analysis to determine the relationships of the different features selected in terms of time frames before seroconversion
	2.9.- DatasetsasinputinPaintomics.R				Creation of tables necessary as inputs for the PaintOmics 3 analysis 
	2.10.-MetabolomicsHeatmapandSUmmaryTable.R			Enrichment analysis of the metabolomics data
	2.11.-VIPs_enrichment_2.R                			Enrichment analysis and Creation of the heatmaps reported in the article
	

#######################
##  Data delivered:  ##
#######################
1.- Folder 	DietaryBiomarkers
	1.-DBprocessed.csv		Dietary Biomarkers Data with the intra-group normalization
	2.-DBRaw.csv			Dietary Biomarkers Data Raw

	Each column represent a time point or a sample individual Id. 
	The rows are:
	Individual.ID: 		ID of the individuals
	Age.in.Months: 		Age in months of the individuals in which that particular sample was evaluated
	Time.to.IA:		It is the timepoint in which the samples were taken relative to the seroconversion point, that is why they are negative numbers, -12 means 12 months before seroconversion and so on.
	From this point on, all rows are the features of the dataset


2.- Folder GeneExpression
	1.-GE_Processed.csv		Gene Expression Data with the intra-group normalization
	2.-GE_Raw.csv			Gene Expression Data Raw

	Each column represent a time point or a sample individual Id. 
	The rows are:
	Sample.Individual.Id:	ID of the individuals samples evaluated
	Individual.ID: 		ID of the individuals
	Age.in.Months: 		Age in months of the individuals in which that particular sample was evaluated
	Time.to.IA:		It is the timepoint in which the samples were taken relative to the seroconversion point, that is why they are negative numbers, -12 means 12 months before seroconversion and so on.
	From this point on, all rows are the features of the dataset

3.- Folder Metabolomics
	Folder Processed
		1.-GCTOF.csv			Gas Chromatography Data with the intra-group normalization
		2.-NegLip.csv			Negative Lipidomics Data with the intra-group normalization
		3.-PosLip.csv			Positive Lipidomics Data with the intra-group normalization

	Folder Raw
		1.-GCTOFRaw.csv			Gas Chromatography Data Raw
		2.-NegLipRaw.csv		Negative Lipidomics Data Raw
		3.-PosLipRaw.csv		Positive Lipidomics Data Raw


	Each column represent a time point or a sample individual Id. 
	The rows are:
	Sample.Individual.Id:	ID of the individuals samples evaluated
	Individual.ID: 		ID of the individuals
	Age.in.Months: 		Age in months of the individuals in which that particular sample was evaluated
	Time.to.IA:		It is the timepoint in which the samples were taken relative to the seroconversion point, that is why they are negative numbers, -12 means 12 months before seroconversion and so on.
	From this point on, all rows are the features of the dataset


4.- AllmetabolitesConverted.txt		Table of functional categories of the all metabolites in study to be able to perform an enrichment analysis
	AllmetabolitesConverted column description:
	ID_Var: 		Id of the variable
	BinBase_name: 		Bin_base name of the metabolite
	Name_abbreviation:	name abbreviated	
	FunCat: 		Functional Category of the metabolite
	FunCatsameName:		Functional Category of the metabolite repeating the name in case they belong to the same category
	Macromolecule:		Macromolecule to which the metabolite belongs


5.- CohortData.csv			Information of the used cohort
	CohortData column description:
	Variable Name:		Definition
	Group.Id:		ID of the matched case/control pair
	Individual.ID: ID of the individuals
	Gender: Male or Female
	First.AAb:		First appearing autoantibody
	Race: White, African American, Unknown/Not Reported/More than one
	Race_ethnicity:  Hispanic,  White Non-Hispanic,  African American, African American,  Missing or Unknown Race and Ethnicity
	Country: USA, Findland, Germany, Sweden
	Case.or.Control: If the individual develops islet of autoimmunity is considered a case
	Model.or.Validation: 	Indicates if the individual was used to create the model (Model) or to validate the model (Validation)

6.-definitivemetabolitesSelectedConvertedmarch9.txt		Table of functional categories of the metabolites selected through VIP_NPLSDA to be able to perform an enrichment analysis
	definitivemetabolitesSelectedConvertedmarch9 column description:
	The same description applies as in AllmetabolitesConverted.txt file


7.-DeltaMatrices_FullArray.RData		List of 4 tables containing the differences in time of the selected variables per individual by period analyzed. It is structured in time frames -12 to -9; -9 to -6, -6 to -3 and -3 to 0. Period 12-9 means the substraction of the values of timepoint -12 minus the timepoint -9 values. The same procedure for eah time point so the results reflects the changes in feature expression as the individuals approach to their seroconversion point. This to be able to create the partial correlation networks. 

8.- LM_globalTargets				Information of the used cohort in a categorical structure, sorted by the sample Individual ID
	LM_globalTargets column description:
	Sample.Individual.ID: 	ID of the individuals samples evaluated
	Individual.ID: 		ID of the individuals
	outcome: 		If the individual develops islet of autoimmunity is considered a case represented with the number 1, if not is a control represented with the number 0
	time: 			It is the timepoint in which the samples were taken
	FirstAAb:		First appearing autoantibody

9.- SelectedFeaturesv3TS			This includes the selected features in each dataset


