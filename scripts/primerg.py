## Primerg (primer3 for gRNAs)
#  Finds NGS primers for multiple gRNA hits within a continuous genomic region
#  Goal:  Find unique NGS primer/ amplicon sequences per targeted site using primer3
#  Input: (Default settings) [1] gRNA list as fasta, [2] continuous genome seq as fasta and 
#                            [3] directory to fasta files 

#  Output: Excel sheet of primers and amplicons per targeted site; 
#          unique ones are marked with "1" in adjacent column

#  Programme outline in default mode: 
#          Process gRNAs in fasta file into a list > List all possibilities of gRNA (+ PAM/PAM-less) >
#          Find position of gRNA (+ PAM/PAM-less)  > Define region for paired-end sequencing (default: 140 bp up/downstream of cleavage site)
#          Design primers with Primer3 within region > Identify unique primers/ amplicons > Filter for unspecific primers (in development)
#          Export excel sheet 

#  Caveats: [1] Amplicons can be derived from copies of targeted genes elsewhere in the genome (not included as input).
#               Problem: filtered reads belong to two genes instead of just the targeted gene.
#               Solution: Check if targeted genes have copies via BLAST. 
#               Concatenate the copies to the continuous genome fasta file supplied to Primerg.

#           [2] Specificity of primers should be evaluated based on reference/ non-reference genome (e.g. via Primer-BLAST)
#               Otherwise, perform two PCRs: first PCR (~1000-1200 bp amplicon) increases specificity of second PCR (~250-280 bp).
#               Note only primers of second PCR is designed in Primerg. 
#               Design first PCR primers by Primer-BLAST (ref. genome) or non-ref genome contigs for increased specificity.

## Requirement(s):
#  Ubuntu
#  Python3 on Linux (primer3 is not available on Windows)
#  Also see 'Import packages' below for required python packages, primarily: primer3, biopython, pandas

## Useful manuals:
#  Arguments for Primer3 command line: http://primer3.org/manual.html

## Variables for users

#  Required input
directory = '/mnt/c/Users/cherw/Desktop/primerg'
gRNA_fasta = "P-loop_minus_set3_gRNA_list.fasta"       #list of gRNAs in 5'-3' direction
genomic_DNA_fasta = "AT1G17615_appended_B3_cluster.fasta"
output_excel = "B3_without_615_set3.xlsx"

#  Optional input

#  Primer3
primer_opt_size = 20                  #recommend 20 bp
primer_min_size = 18                  #recommend 18 bp
primer_max_size = 30                  #recommend 30 bp
primer_opt_tm = 65                    #recommend 65
primer_min_tm = 50                    #recommend 50
primer_max_tm = 70                    #recommend 70
primer_min_GC = 35                    #recommend 40%      *****************************************    
primer_max_GC = 60                    #recommend 60%
primer_product_size_range = [250, 280] #recommend [250, 280] 
primer_dna_conc = 50                  #in nm
primer_dNTP_conc = 0.5                #in mM
primer_pair_max_diff_TM = 5           #recommend 5
primer_salt_divalent = 1.5            #mM of divalent cations e.g. Mg2+
primer_num_return = 100               #maximum number of primers to return per primer3 run

#  CRISPR-Cas parameters
pam = 'NG'                           #e.g. "NGG", let pam be False if pamless ****************************************
gRNA_len = 20
cleavage_pos = 3                      #nth base from 3' end of gRNA where Cas9 cleaves genomic template,
                                      #default = 3 where NNNNNNNNNNNNNNNNN|NNN 

#  Primerg
simplified_output = False              #if True, output file will only include unique primers/ amplicons. If False, output file will consist of all primers/ amplicons.


#######################################################################################################################
## Source code
#  Do not edit unnecessarily

## Import packages
import os
import primer3
from Bio import SeqIO
from Bio.Seq import Seq
import re
import pandas as pd

## Reading parameters
os.chdir(directory)
template = str(SeqIO.read(genomic_DNA_fasta, "fasta").seq).upper()    #template read as string


## Read gRNA fasta file into a list
def gRNA_listing(gRNA_fasta):
    gRNA_list = []
    try:
        gRNA_fasta_list = SeqIO.parse(gRNA_fasta, "fasta")
        for record in gRNA_fasta_list:
            gRNA_list += [str(record.seq).upper()]      
    except ValueError:                #to catch fasta files with single gRNA (use SeqIO.read instead of SeqIO.parse)
        gRNA_fasta_list = SeqIO.read(gRNA_fasta, "fasta")
        gRNA_list += [str(gRNA_fasta_list.seq.upper())] 
    return gRNA_list


## List all possibilities of gRNA + pam
def gRNA_pam_listing(gRNA_list, pam, template):
    gRNA_pam_list_F = []
    gRNA_pam_list_R = []
    gRNA_pam_list = []    
    
    #for CRISPR with pam
    if pam:
        #find all possibilities of pam    
        pam_list = [pam.replace('N', 'A'), pam.replace('N', 'T'), pam.replace('N', 'C'), pam.replace('N', 'G')]
        
        #generate list of gRNA + pam (forward)
        temp_list_F = [a + s for s in pam_list for a in gRNA_list] 
        
        #generate list of gRNA + pam (reverse complement)
        temp_list_R = [str(Seq(s).reverse_complement()).upper() + str(Seq(a).reverse_complement()).upper() for s in pam_list for a in gRNA_list]
        
        #create new list of gRNA + pam that is found in template
        for gRNA_pam in temp_list_F:
            if gRNA_pam in template:
                gRNA_pam_list_F += [ gRNA_pam ]
        for gRNA_pam in temp_list_R:
            if gRNA_pam in template:
                #converts gRNA + pam into reverse complement before adding to list
                gRNA_pam_list_R += [str(Seq(gRNA_pam).reverse_complement())]
        #output list will consist of [ [forward gRNA + pam], [reverse complement gRNA + pam]]
        gRNA_pam_list = [gRNA_pam_list_F, gRNA_pam_list_R]
        return gRNA_pam_list
    
    #for pamless CRISPR
    elif pam == False: 
        #generate list of gRNA (forward)
        temp_list_F = [a for a in gRNA_list] 
        
        #generate list of gRNA (reverse complement)
        temp_list_R = [str(Seq(a).reverse_complement()).upper() for a in gRNA_list]
        
        #create new list of gRNA that is found in template
        for gRNA_pam in temp_list_F:
            if gRNA_pam in template:
                gRNA_pam_list_F += [ gRNA_pam ]
        for gRNA_pam in temp_list_R:
            if gRNA_pam in template:
                #converts gRNA + pam into reverse complement before adding to list
                gRNA_pam_list_R += [str(Seq(gRNA_pam).reverse_complement())]
        #output list will consist of [ [forward gRNA + pam], [reverse complement gRNA + pam]]
        gRNA_pam_list = [gRNA_pam_list_F, gRNA_pam_list_R]
        return gRNA_pam_list        

## Primer design

#  Create error class
class SequenceError(Exception):
    pass

#  Find gRNA in genomic DNA: returns list of positions of each gRNAs; only one position per gRNA will be returned
def gRNA_finder(gRNA_pam_list, template):
    temp = []
    lst_F = []
    lst_RC = []
    lst_grand = []
    
    #note: data structure of gRNA_pam_list as [ [forward gRNA + pam], [reverse complement gRNA + pam]]
    
    #find position of forward gRNA cleavage site
    for gRNA_pam in gRNA_pam_list[0]:
        #find overlapping matches for forward gRNA AND keep list of Cas9 cleavage site position
        temp  += [ (m.start() + len(gRNA_pam) - len(pam) - cleavage_pos) for m in re.finditer('(?=' + gRNA_pam + ')', template) ] 
        lst_F += [(gRNA_pam, temp)]
        temp = []
        
    #find position of reverse complement gRNA cleavage site
    for gRNA_pam in gRNA_pam_list[1]:
        #find overlapping matches with reverse comp gRNA AND keep list Cas9 cleavage site position
        temp += [ (m.start() + len(pam) + cleavage_pos) for m in re.finditer('(?=' + str(Seq(gRNA_pam).reverse_complement()) + ')', template) ] 
        lst_RC += [(gRNA_pam, temp)]
        temp = []
        
    #structure of lst_grand is e.g. [('GGCGTTGACAGATGAGGGGCAGG', [403, 501]), ('AATGCTGGATTTTCTGCCTGTGG', [643])]
    lst_grand += lst_F + lst_RC
    return lst_grand

    
#  Design primers      
                
def primer_design(gRNA_pos, template):
    #initiate dataframe
    df = pd.DataFrame(columns = ['gRNA_index', 'region_index', 'combined_index', 'gRNA_pam', 'region', 'primer_F', 'unique_primer_F', 'primer_R', 'unique_primer_R', 'amplicon_F', 'unique_amplicon_F', 'amplicon_R', 'unique_amplicon_R'])                             
    
    #find targeted regions where primers are designed
    half_len = round(primer_product_size_range[1]/2) #half the length of amplicon
    gRNA_index = 0
    region_index = 1
    
    for unit in gRNA_pos:
        gRNA_seq = unit[0]
        pos_lst = unit[1]
        gRNA_index += 1
        #iterate through the position of cleavage sites
        for pos in range(len(pos_lst)):
            
            #reset regions for this cleavage site 
            regions = []
            
            #initiate template seqs per cleavage site
            regions += [template[pos_lst[pos] - half_len: pos_lst[pos] + half_len]]
            
            #design primers for current position
            for i in regions:
                design_primer = primer3.bindings.designPrimers(
                    {
                        'SEQUENCE_TEMPLATE': i
                    },
                    {
                        'PRIMER_OPT_SIZE': primer_opt_size,
                        'PRIMER_MIN_SIZE': primer_min_size,
                        'PRIMER_MAX_SIZE': primer_max_size,
                        'PRIMER_OPT_TM': primer_opt_tm,
                        'PRIMER_MIN_TM': primer_min_tm,
                        'PRIMER_MAX_TM': primer_max_tm,
                        'PRIMER_MIN_GC': primer_min_GC,
                        'PRIMER_MAX_GC': primer_max_GC,
        
                        'PRIMER_DNA_CONC': primer_dna_conc,
                        'PRIMER_PRODUCT_SIZE_RANGE': primer_product_size_range,
                        'PRIMER_DNTP_CONC': primer_dNTP_conc,
                        'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                        'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                        'PRIMER_NUM_RETURN': primer_num_return,
                        
                    })
                
                #retrieve all primers from design 
                for j in range(len(design_primer)):
                    left = None
                    right = None
                    try:
                        left = design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE']
                        right = design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE']
                        
                        #dataframe format: ['gRNA', 'region', 'primer_F', 'primer_R', 'amplicon_F', 'amplicon_R', 'set')
                        #index will refer to unique gRNA + region combination
                        
                        #find reverse primer position at 5' end (back)
                        primer_R_pos = template.find(str(Seq(right).reverse_complement())) + len(right)
                        
                        df2 = pd.DataFrame([ [gRNA_index, region_index, f"{gRNA_index}_{region_index}", gRNA_seq, i, left, '', right, '', 
                                              template[template.find(left): template.find(left) + half_len], '', #forward amplicon
                                              template[primer_R_pos - half_len : primer_R_pos], ''],],  #reverse amplicon in 5'-3'
                                           columns = ['gRNA_index', 'region_index', 'combined_index', 'gRNA_pam', 'region', 'primer_F', 'unique_primer_F', 'primer_R', 'unique_primer_R', 'amplicon_F', 'unique_amplicon_F', 'amplicon_R', 'unique_amplicon_R'])
                        df = df.append(df2)    
                    
                    except KeyError:
                        df2 = pd.DataFrame([ [gRNA_index, region_index, f"{gRNA_index}_{region_index}", gRNA_seq, i, 'NA', 'NA', 'NA', 'NA', 'NA' , 'NA', 'NA', 'NA'],],
                                           columns = ['gRNA_index', 'region_index', 'combined_index', 'gRNA_pam', 'region', 'primer_F', 'unique_primer_F', 'primer_R', 'unique_primer_R', 'amplicon_F', 'unique_amplicon_F', 'amplicon_R', 'unique_amplicon_R'])
                        df = df.append(df2) 
                        break
                #increase count for region_index
                region_index += 1                    
    return df

##  Find unique primers and amplicons between groups
#   Count of primers should be non-overlapping (because overlapping primer sequnces within ~30 bp range can be considered as the same primers)
#   Count of amplicons should be overlapping (because generating precise unique amplicon is important to filter reads to correct targeted sites accordingly)

def occurrences(substring, string):
    #count number of times a substring occurs in a string 
    #overlapping substrings are counted
    return len(re.findall('(?={0})'.format(re.escape(substring)), string))  

def rev_complement(string):
    #reverse complement a string and return a string
    return str(Seq(string).reverse_complement())

def find_unique(df):
    df_temp = None
    
    #identify unique primers/ amplicons between groups
    for index, row in df.iterrows():
        #create separate dataframe of other groups        
        current_index = row['combined_index']
        df_temp = df[df.combined_index != current_index]
        
        
        ##check for unique amplicon F 
        if occurrences(row['amplicon_F'], template) == 1 and occurrences(row['amplicon_F'], rev_complement(template)) == 0:
            row['unique_amplicon_F'] = 1
            
            #check if primer F is not found in primer list from other targeted sites AND unique in genome
            if row['primer_F'] not in (list(df_temp['primer_F']) + list(df_temp['primer_R'])):
                if occurrences(row['primer_F'], template) == 1 and occurrences(row['primer_F'], rev_complement(template)) == 0:
                    row['unique_primer_F'] = 1
                elif occurrences(row['primer_F'], template) > 1 or occurrences(row['primer_F'], rev_complement(template)) > 0:
                    row['unique_primer_F'] = 0
            else:
                row['unique_primer_F'] = 0
            
        elif occurrences(row['amplicon_F'], template) > 1 or occurrences(row['amplicon_F'], rev_complement(template)) > 0:
            row['unique_amplicon_F'] = 0
            row['unique_primer_F'] = 0
            
        elif occurrences(row['amplicon_F'], template) == 0:
            row['unique_amplicon_F'] = "Error, amplicon is not found in genomic template"
            row['unique_primer_F'] = "Error, primer is not found in genomic template"
            
            
        ##check for unique amplicon R
        if occurrences(row['amplicon_R'], template) == 1 and occurrences(row['amplicon_R'], rev_complement(template)) == 0:
            row['unique_amplicon_R'] = 1
            
            #check if primer R is not found in primer list from other targeted sites
            if row['primer_R'] not in (list(df_temp['primer_F']) + list(df_temp['primer_R'])):
                if occurrences(row['primer_R'], rev_complement(template)) == 1 and occurrences(row['primer_R'], template) == 0:
                    row['unique_primer_R'] = 1
                elif occurrences(row['primer_R'], rev_complement(template)) > 1 or occurrences(row['primer_R'], template) > 0:
                    row['unique_primer_R'] = 0
            else:
                row['unique_primer_R'] = 0
        
        elif occurrences(row['amplicon_R'], template) > 1 or occurrences(row['amplicon_R'], rev_complement(template)) > 0:
            row['unique_amplicon_R'] = 0
            row['unique_primer_R'] = 0
            
        elif occurrences(row['amplicon_R'], template) == 0:
            row['unique_amplicon_R'] = "Error, amplicon is not found in genomic template"
            row['unique_primer_R'] = "Error, primer is not found in genomic template"
            
    return df

def export(df_unique):
    #drop rows without any unique primers/ amplicons
    if simplified_output == True:
        #reset index to get index from 0 to end to give each row a unique index
        df_unique.reset_index(drop=True, inplace=True)
        
        #iterrate through the rows and drop rows accordingly
        indexes_to_drop = []
        for index, row in df_unique.iterrows():
            if row['unique_primer_F'] == row['unique_primer_R'] == row['unique_amplicon_F'] == row['unique_amplicon_R'] == 0:
                indexes_to_drop += [index]
        df_unique.drop(df_unique.index[indexes_to_drop], inplace=True)
    
    #export dataframe as excel sheet
    df_unique.to_excel(output_excel, index = False)
    pass
     

def execute_all():
    gRNA_list = gRNA_listing(gRNA_fasta)
    print("gRNA list obtained.")
    
    gRNA_pam_list = gRNA_pam_listing(gRNA_list, pam, template)
    if pam:
        print("gRNA + PAM list obtained.")
    elif pam == False:
        print("PAM-less gRNA mode active.")
        
    gRNA_pos = gRNA_finder(gRNA_pam_list, template)
    if pam:
        print("All gRNA + PAM positions are found in template.")
    elif pam == False:
        print("All PAM-less gRNA positions are found in template.")
    
    print("")
    print("Designing primers...")
    df = primer_design(gRNA_pos, template)
    print("Primers designed and amplicons are retrieved.")
    
    print("Finding unique primers and amplicons...")
    df_unique = find_unique(df)
    print("Done.")
    
    print("")
    export(df_unique)
    print(f"Output is placed at {directory} as {output_excel}.")

##  Run programme
execute_all()

#Feedback
#2. Filter to keep only rows with 1(s)
#3. Primer BLAST!
#4. template.count for amplicon should consider overlapping but not for primer
#5. find signature motif for read sorting: only for unique amplicon = 1 but unique primer = 0




