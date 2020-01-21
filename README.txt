Program: Normalizing Counts 
Created by M. Joseph Tomlinson

Note: Code was Developed in the Abasht Laboratory at the University of Delaware 
under the supervision of Dr. Behnam Abasht website: "http://canr.udel.edu/faculty/behnam-abasht/"

Program normalizes the counts of allele specific expression (ASE) based on the number of biallelic variants for a gene. This tries
to take into account larger genes which may show more ASE, but when overall number of biallelic variants for
that gene are accounted for, it may show a lower percentage of ASE.

Equation on a Per Gene Basis (significant ASE counts/total counts) = Percentage ASE per gene

Program utilizes information from "Ensemble_Gene_IDs" count files created by the VEP_Results_Parser_Merger program and also
uses an annotated variant file created by Ensembl's VEP program (https://useast.ensembl.org/info/docs/tools/vep/index.html)

Also reported by the program is the exon and intron counts for only the significant variants for each gene identified.


Program Requirements:
python 3.6


How to Run Program:
Update the parameter file named (Normalizing_Counts_Parameter_File.txt) with the required inputs and double click on program to run


Required Input Files:
Importnat note: The required input files were created by Ensemble's VEP and VEP_Results_Parser_Merger Program
prior utilized in the analysis of data.

Three Files:

1) Testable_Counts_File
- File consists of all testable genes with overall counts reported by the VPRM program. It is recommended to use
the Ensemble gene IDs counts file.

2) Significant_Counts_File
- File consists of genes that showed statistically significant ASE variants reported by VRPM Program. It is recommended to use
the Ensemble gene IDs counts file.

3) Annotation File (Sig Variants Only)- VEP Output File
- Files of all variants and their corresponding biological impact information as annotated by VEP.

Program Output:

1) normalized_counts_values.txt
- File consists of all the normalized counts and overall counts of the exons and introns for significant variants

2) summary_report_file.txt
- A summary report file of basic information about the data run (files used) and overall basic statistics of the run

3) not_sig_genes.txt
- All genes not found in the significant genes file, meaning genes showed no sign of ASE. 




