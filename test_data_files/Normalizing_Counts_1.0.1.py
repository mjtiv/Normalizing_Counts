#!/usr/bin/env python3.6


"""

Program: Normalizing Counts
Created by M. Joseph Tomlinson

Note: Code was Developed in the Abasht Laboratory at the University of Delaware 
under the supervision of Dr. Behnam Abasht website: "http://canr.udel.edu/faculty/behnam-abasht/"

Program normalizes the counts of ASE based on the number of biallelic variants for a gene. This tries
to take into account larger genes may show more ASE, but when overall number of biallelic variants for
that gene are accounted for, it may show a lower percentage of ASE.

Equation on a Per Gene Basis (significant ASE counts/total counts) = Percentage ASE per gene

Program utilizes information from "Ensemble_Gene_IDs" count files created by the VEP_Results_Parser_Merger program.

Also reported by the program is the exon and intron counts for only the significant variants for each gene identified. 


"""




def parsing_input_parameter_file():
    
    """ 

    Parses through the paramter file to return a dictionary
    of the input parameters
    
    :param none: automatically opens file

    :return dictionary: dictionary of file names for analysis
    
    """

    # Input parameter file being parsed
    input_file = open('Normalizing_Counts_Parameter_File.txt', 'r')

    # Summary Report file for record keeping
    summary_report_file = open('summary_report_file.txt', 'w')

    print ("Parsed Lines from User")
    
    for line in input_file:
        if line.startswith("--"):
            line=line.rstrip('\n')
            # just printing user inpt for validation
            print (line)
            parsed_parameters=line.split("--")

            for x in range(1, len(parsed_parameters)):
                inputs = parsed_parameters[x].split(" ")
                
                if inputs[0] == "Testable_Counts_File":
                    testable_counts_file = inputs[1] 

                elif inputs[0] == "Significant_Counts_File":
                    significant_counts_file = inputs[1]

                elif inputs[0] == "Sig_Var_VEP_Annot_File":
                    sig_var_vep_annot_file = inputs[1]

                else:
                    print ("Incorrect file submitted")

    # Skip anything else
    else:
        pass

    # Printing to Command Line for Troubleshooting
    print ("Name of Testable Counts File: ", testable_counts_file)
    print ("Name of Significant Var. Counts File: ", significant_counts_file)
    print ("Name of Significant Var. Annotation File: ", sig_var_vep_annot_file)

    # Writing to an actual file
    summary_report_file.write("Running Normalization of Data Program")
    summary_report_file.write("\n")
    summary_report_file.write("Name of Testable Counts File: " + str(testable_counts_file) + "\n")
    summary_report_file.write("Name of Significant Counts File: " + str(significant_counts_file) + "\n")
    summary_report_file.write("Name of Significant Var. Annotation File: " + str(sig_var_vep_annot_file) + "\n")
    summary_report_file.write("\n")

    summary_report_file.close()
    input_file.close()
    
    return{'testable_counts_file':testable_counts_file,
           'significant_counts_file':significant_counts_file,
           'sig_var_vep_annot_file':sig_var_vep_annot_file}

    
def parse_file_create_dict(file_name):
    
    """

    Opens up file and parses through data creating
    a dictionary of results
    
    :param file_name: name of file being parsed
    :return results_dict: dictionary of parsed results
    
    """

    # create empty dictionary
    results_dict = {}
    countable_values = 0
    
    input_file = open(file_name, 'r')

    for line in input_file:
        # Skips the no_ID during normalization
        if line.startswith(('Ensemble', 'no_ID')):
            continue

        else:
            line = line.rstrip('\n')
            split_data = line.split("\t")
        
            results_dict.update({split_data[0]: split_data[1]})

            # Tally number of values (exclude no hits)
            countable_values += 1 

    input_file.close()

    return (results_dict)


def parse_vep_results(vep_results_file):

    '''

    Parses through the VEP output to identify the number of
    unique exon locations of variants in the file, which will
    later be merged back with the normalized data, to help prioritize
    ASE locations

    : Param vep_results_file: 
    
    '''
    
    input_file = open (vep_results_file, 'r')

    vep_location_dict = {}

    for line in input_file:

        if line.startswith("#"):
            continue
        else:
            line = line.split('\t')

            # Get ensemble gene ID
            ensemble_id = line[5]

            # Get exon value
            exon_location = line[9]
            intron_location = line[10]

            if exon_location == "-" and intron_location == "-":
                continue

            if ensemble_id in vep_location_dict:
                # taking into account exons and intron designations
                if intron_location == "-":
                    vep_location_dict[ensemble_id]['exon'].append(exon_location)

                else:
                    vep_location_dict[ensemble_id]['intron'].append(intron_location)

            # Creating a new value in the dictionary
            else:
                
                if intron_location == "-":
                    vep_location_dict.update({ensemble_id: {'exon': [exon_location], 'intron': []}})

                # Working Here  
                else:
                    vep_location_dict.update({ensemble_id: {'exon': [], 'intron': [intron_location]}})
              
    input_file.close()

    return(vep_location_dict)


def normalize_values(testable_dict, significant_hits_dict, sig_var_annot_dict):

    output_results_file = open('normalized_counts_values.txt', 'w')
    summary_report_file = open('summary_report_file.txt', 'a')
    non_sig_values_list_output = open('not_sig_genes.txt', 'w')

    output_results_file.write("Ensemble_ID\tTestable_Count\tSig_Counts\tNormailized_Value(%)\tTotal_Sites" +
                              "\tExon_Sites\tIntron_Sites\n")

    # Setup Counters
    tested_values_counter = 0
    match_counter = 0
    no_matches_counter = 0
    non_matches_list = []

    for key in testable_dict:
        tested_values_counter += 1

        try:
            key_value = key
            testable_count = int(testable_dict[key])
            significant_count = int(significant_hits_dict[key])
            normalized_value_percent = str(round((significant_count / testable_count * 100),2))

            # Getting splicing location information
            try: 
                exon_sites = set(sig_var_annot_dict[key_value]['exon'])
                intron_sites = set(sig_var_annot_dict[key_value]['intron'])
                total_exon_intron_sites= len(exon_sites) + len(intron_sites)

                output_results_file.write(str(key) + "\t" + str(testable_count) + "\t"
                                          + str(significant_count) + "\t"
                                          + normalized_value_percent + "\t"
                                          + str(total_exon_intron_sites) + "\t"
                                          + str(len(exon_sites)) + "\t"
                                          + str(len(intron_sites)) + "\n")

            # No Exon or Intron Sites Reported
            except KeyError:
                exon_sites = "-"
                intron_sites = "-"
                total_splice_sites= "-"
            
                output_results_file.write(str(key) + "\t" + str(testable_count) + "\t"
                                          + str(significant_count) + "\t"
                                          + normalized_value_percent + "\t"
                                          + str(total_exon_intron_sites) + "\t"
                                          + str(exon_sites) + "\t"
                                          + str(intron_sites) + "\n")
            match_counter += 1

        except KeyError:
            non_sig_values_list_output.write(str(key) + "\n")
            no_matches_counter += 1
            non_matches_list.append(key)
            continue

    # Write summary stats to summary report file
    summary_report_file.write("The total number of values in testable file are: "
                              + str(tested_values_counter) + "\n")
    summary_report_file.write("The total number of matching values between files are: "
                              + str(match_counter) + "\n")
    summary_report_file.write("The total number of non-matching values between files are: "
                              + str(no_matches_counter)+ "\n" )

    # Close the files
    non_sig_values_list_output.close()
    output_results_file.close()
    summary_report_file.close()
    

    return()

def main():

    #########################################################
    print ("Starting Run of the Program")
    
    # Opens the parameter file to get file name
    print("Opening the parameter file")
    parameter_stuff = parsing_input_parameter_file()

    # Retrieving the parameters from the parameter file parsing output
    testable_input_file = parameter_stuff['testable_counts_file']
    significant_input_file = parameter_stuff['significant_counts_file']
    sig_var_vep_annot_file = parameter_stuff['sig_var_vep_annot_file']
    
    print ("Program retrieving data from input files")

    testable_dict = parse_file_create_dict(testable_input_file)

    significant_hits_dict = parse_file_create_dict(significant_input_file)

    sig_var_annot_dict = parse_vep_results(sig_var_vep_annot_file)

    print ("All data retrieved succesfully from input files")

    print ("Normalizing Values")
    normalize_values(testable_dict, significant_hits_dict, sig_var_annot_dict)

    
    print ("Program Done Running")
main()


# Version Control
# Normalizing_Counts_1.0.1.py
# - Annotated code better and changed the "splice_sites" variable name to "total_exon_intron_sites," overall naming
# of variable was poorly worded initially
# 
# Normalizing_Counts_1.0.0.py
# -Fixed Bug where no exons or splice sites for variant calls were being discarded


