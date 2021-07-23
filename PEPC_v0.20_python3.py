# -*- coding: utf-8 -*-
"""
Created on Wed Dec 05 15:04:04 2012

@author: Jeff Bowman, bowmanjs@gmail.com

v1.15 adds secondary structure using predictions from psipred.  mutations only
occur in the selected secondary structure elements, or in all elements if 
secondary structure is not specified.  It might be desirable in a future version
to allow secondary structures not selected to mutate in a random fashion.  However
since there is no way to discriminate between regions of high and low conservation
this may overrepresent the number of random mutations taking place

v1.16 adds uses hmmscan against Pfam-A to detect conserved domains.  These
domains are masked from mutation.

v0.17 resets the version numbers to reflect development stage.  This version
also allows PEPC to use multiple parameters and multiple secondary structure
elements to select survivors.

v0.18 allows only one mutation per iteration.  Iterations with rare positive
mutations were capturing too many neutral mutations.  Also elimate some ugly
nested dictionaries, now using paired keys.  This version also corrects a
serious error where non-beneficial mutations were surviving due to incorrect
reassignment of the aa_seq object (see http://stackoverflow.com/questions/19341365/setting-two-arrays-equal).

v0.19 impliments Sarah Purkey's suggestion for improving residues to mutate
and candidate residues.

v0.20 uses blastp to estimate a rate parameter for each position, instead of
using the binary rate provided by assuming a conserved domain never mutates. 

v0.21 removes the influence of residue composition on inital rate variation parameter.  The
idea here is to isolate variation based on position alone, since variation by
residue composition is accounted for later.  Note that when viewing a histogram
of the posit_var.txt file you probably won't be able to distinguish conserved
domains (you could in v20), because the conserved nature of the domain is also
dependent on amino acid composition, which is not reflected at that step.

v0.21 is DEFUNCT - Was radically underestimating the alpha (shape) parameter
of rate variation.  Reverting to v20 which produces realistic alpha parameter.

To do: write output directly to file at each mutation, no need to save to nested list
currently each generation results in multiple outputs anyway, one for each mutation made

"""

#### user setable parameters ####

db = 'Pfam-A.fasta' ## database for blastp search

#### end user setable parameters ####

import random
from Bio.SubsMat import MatrixInfo
import math
from Bio.SeqUtils import ProtParam as PP
from Bio import SeqIO
import sys
import re
import subprocess
import os
import string

# print info and copyright

print('This is the Protein Evolution Parameter Calculator (PEPC)\n\
Copyright Jeff Bowman, 2012')

# provide help as an option

help = 'Run PEPC as\n\
python PEPC_vX.py [options] [fasta] [starts] [events] [dice]\n\
\n\
where:\n\
X = the version number\n\
options = a combination of secondary structure elements coil (C), beta sheets\n\
(E), or alpha helices (H) and parameters (flex, gravy, iso, arom, ai) and\n\
direction (up, down), e.g. C:flex=up,gravy=down.\n\
fasta = a protein fasta in the working directory containing a single\n\
sequence, omit the .fasta extension.  The file name can NOT have a colon in it.\n\
Sorry.\n\
starts = the number of replicate runs\n\
events = the number of mutations you would like to attempt to see an\n\
improvement in the parameter space.  Setting to 0 will calculate the\n\
protein parameters by secondary structure.\n\
dice = the number of sides to the dice rolled to determine if a mutation that\n\
does not improve fitness should survive.  0 turns off dice, meaning only\n\
mutations that improve fitness will be kept.\n\
\n\
PEPC outputs four sets of files.  For each start a new fasta is created\n\
containing the mutated aa sequence.  For each program run an output file\n\
containing the protein parameters for each mutation is created.  This output\n\
file is in the order: start, event (0 = the unmutated sequence), isoelectric\n\
point, flexibility, aromaticity, gravy, aliphatic index.  The tally.txt file\n\
contains the mutations by position.  A log file is also generated for each\n\
program run.  The line is printed for each event in the order start, event,\n\
number of mutations, number of random mutations, position of possible mutation,\n\
current residue at that position, candidate residue, success of mutation, and\n\
resulting change to the indicated parameters.\n\
\n\
Protein parameter isoelectric point, flexibility, aromaticity, and gravy are\n\
calculated using the ProtParam.ProteinAnalysis module in Biopython.\n\
Aliphatic index is calculated using the following equation:\n\
100*[A] + 2.9*100*[V] + 3.9*100*[I]+100*[L] according to Akai, 1980.\n\
\n\
In addition to Biopython PEPC requires the use of psipred.  You will need to \n\
download and install this dependency separately, along with your preferred\n\
protein database (recommend Pfam-A.fasta).  Specifically PEPC makes use\n\
of the runpsipredplus script, so make sure that this script is properly\n\
configured and available in your path.  I recommend running the PSI-Blast\n\
against Pfam-A.fasta, after conversion this file to a blast database.  Modify\n\
the blast line of the runpsipredplus file accordingly.  You may also wish to\n\
specify the number of cpus for blast to use with the num_threads parameter on\n\
this line.  Refer to the blast documentation at\n\
http://www.ncbi.nlm.nih.gov/books/NBK1763/ for more details.\n\
\n\
To determine site-specific mutation rates PEPC relies on blastp.  It is\n\
recommended that you run blastp against the Pfam-A database used by psipred, \n\
but you need to tell PEPC where to look for this database by setting the\n\
\'db\' parameter at the top of the PEPC python file.\n\
\n\
Full command example:\n\
python PEPC_v0.17.py C:flex=down,gravy=up E:arom=up fasta_file_basename 3 1000 10000'

#### TEST SWITCH ####
test = False ## Set to false if you want command line arguments used
#####################

test_list = ['null path', 'C:flex=up', 'H:flex=up,iso=down', 'beta_galactosidase', '1', '10', '0'] ## a list used for testing or interactive work

if '-h' in sys.argv:
    print(help)
    sys.exit()

## create new list of options and purge options from sys.argv

required = []
direction = [] ##  [C,E,H]:param:direction with capitalization as shown
select = '' ## used for labeling output files

#### parse command line or test list and define some variables ####

if test == False:
    for value in sys.argv:
        s = str(value)
        if ':' in s:
            direction.append(s)
            select = select + '_' + s
        else:
            required.append(value)

elif test == True:
    for value in test_list:
        s = str(value)
        if ':' in s:
            direction.append(s)
            select = select + '_' + s
        else:
            required.append(value)

select = str.replace(select, ':', '_')
select = str.replace(select, ',', '_')

## set variables based on command line entries

exp_seq = required[1]
starts = int(required[2])
events = int(required[3])

if int(required[4]) == 0:
    r_dice = -1 ## -1 not in population, so won't be selected
    dice = range(0,int(1),1) ## you have to give it something, or returns sampling error later
else:
    dice = range(0,int(required[4]),1)
    r_dice = random.sample(dice, 1)[0]

## parse the direction list from the command line
    
director = {} ## key = SS, value = {param:direction, ...}

for each in direction:
    each = each.split(':')
    ss = each[0]
    param_direct = each[1]
    param_direct = param_direct.split(',')
    
    temp_dict = {}
    
    for pd in param_direct:
        pd = pd.split('=')
        param = pd[0]
        direct = pd[1]
        temp_dict[param] = direct
     
    director[ss] = temp_dict
    
## read in fasta and matrix

exp_record = SeqIO.read(exp_seq+'.fasta','fasta')
aa_seq_str_i = str(exp_record.seq)

sub_matrix = MatrixInfo.blosum80

## must clean B, X and Z from matrix

for entry in list(sub_matrix):
    if entry[0] == 'B':
        del sub_matrix[entry]
    elif entry[1] == 'B':
        del sub_matrix[entry]
    elif entry[0] == 'Z':
        del sub_matrix[entry]
    elif entry[1] == 'Z':
        del sub_matrix[entry]
    elif entry[0] == 'X':
        del sub_matrix[entry]
    elif entry[1] == 'X':
        del sub_matrix[entry]
        
#### function definitions ####
        
## mutation matrix function
        
def build_mut_matrix(seq, aa_freq, sub_matrix, posit_var):
    
    probs = {} ## key = (position, residue) value = likelihood of mutation occurring

    for i,aa in enumerate(seq):
        try:
            i_posit_var = posit_var[i]
        except KeyError:
            i_posit_var = 1
            
        for residue in aa_freq.keys():
            if aa != residue:
                freq = aa_freq[aa] * aa_freq[residue]
                try:
                    prob = freq * math.pow(10, sub_matrix[aa, residue]) * i_posit_var
                except KeyError:
                    prob = freq * math.pow(10, sub_matrix[residue, aa]) * i_posit_var
            else:
                prob = 0
            probs[i, residue] = prob
    
    all_sums = []        
    sums = 0
    sorted_keys = sorted(probs.keys())
    
    for key in sorted_keys:
        sums = sums + probs[key]
        probs[key] = sums
        all_sums.append(sums)
                
    min_sum = min(all_sums)
    
    for key in probs.keys():
        try:
            probs[key] = int(probs[key] / min_sum)
        except ZeroDivisionError:
            continue
            
    return probs, sorted_keys, probs[sorted_keys[-1]]
    
## function to select a mutation candidate
    
def get_cand(sorted_keys, probs, seq, sums):
            
    r = random.randint(0, int(sums))    
    for key in sorted_keys:
        if probs[key] >= r:
            break
        
    return key[0], seq[key[0]], key[1]
        
## diagnostic function, remove when complete

def diagnostic(iter1, iter2):
    p = 0
    for i,j in enumerate(iter1):
        if j != iter2[i]:
            p = p + 1
    print(p)
    return p
    
## define function to calculate protein parameters

def calculate_pp(sequence):
    seq = PP.ProteinAnalysis(sequence)
    iso = PP.ProteinAnalysis.isoelectric_point(seq)
    flex_list = PP.ProteinAnalysis.flexibility(seq)
    
    if len(flex_list) == 0:
        flex = 'NA'
    else:
        flex = sum(flex_list)/len(flex_list)
        
    arom = PP.ProteinAnalysis.aromaticity(seq)
    gravy = PP.ProteinAnalysis.gravy(seq)
    aa = PP.ProteinAnalysis.get_amino_acids_percent(seq)
    ai = (100 * aa['A']) + (2.9 * 100 * aa['V']) + (3.9 * (100 * aa['I']+100 * aa['L']))  
    
    return flex, gravy, iso, arom, ai
    
## define function to handle writing nested list as csv
    
def nested_csv(handle, list_name):
    import csv
    nested_list = csv.writer(handle, dialect='excel-tab', lineterminator='\n')
    nested_list.writerows(list_name)  
    
#### set some initial variabiles ####
    
pos_w_mut = [] ## number of mutations by position
param_output = [] ## parameter calculations will go here
coil_param_output = []
beta_param_output = []
helix_param_output = []

## define amino acid frequencies for log likelihood calculations

aa_freq = {}
aa_freq['A'] = 0.075
aa_freq['R'] = 0.052
aa_freq['N'] = 0.046
aa_freq['D'] = 0.052
aa_freq['C'] = 0.018
aa_freq['E'] = 0.041
aa_freq['Q'] = 0.063
aa_freq['G'] = 0.071
aa_freq['H'] = 0.022
aa_freq['I'] = 0.055
aa_freq['L'] = 0.091
aa_freq['K'] = 0.058
aa_freq['M'] = 0.028
aa_freq['F'] = 0.039
aa_freq['P'] = 0.051
aa_freq['S'] = 0.074
aa_freq['T'] = 0.060
aa_freq['W'] = 0.013
aa_freq['Y'] = 0.033
aa_freq['V'] = 0.065

## clean sequence string

aa_seq_str_i = aa_seq_str_i.replace('\n','')
aa_seq_str_i = aa_seq_str_i.replace('X','')

temp_fasta = open(exp_seq+'.clean.fasta', 'w')
print( '>temp\n '+ aa_seq_str_i, file = temp_fasta)
temp_fasta.close()


#### calculate rate parameter ####

## run blast

if exp_seq + '.blast' not in os.listdir('.'):
    print('running blastp to determine rate variation in query sequence')
    os.system('/home/lpiszkin/PEPC/bin/sh/blastp -num_alignments 1000 -evalue 1e-5 -query \"'+exp_seq+'.clean.fasta\" -db ' + db + ' -outfmt 0 -out \"' + exp_seq + '.blast\"')

## parse blast output

posit_var = {} ## key = position with start = 0, value = [number of alignments, number of mismatches]

with open(exp_seq + '.blast', 'r') as blast_in:
    record = False
    for line in blast_in:
        if line.startswith('>'):
            record = True
        
        if record == True:
            if line.startswith('Query'):
                line = line.split()
                qstart = int(line[1])
                qstring = line[2]
            
            elif line.startswith('Sbjct'):
                line = line.split()
                sstart = line[1]
                sstring = line[2]
        
                for i,r in enumerate(qstring):
                    try:
                        temp = posit_var[i + qstart - 1]
                    except KeyError:
                        temp = [0, 0]
                    temp[0] = temp[0] + 1
                    if r != '-':
                        if r != sstring[i]:
                            temp[1] = temp[1] + 1
                    posit_var[i + qstart - 1] = temp

## convert to fraction of alignments that differ from query
                   
for key in posit_var.keys():
    temp = posit_var[key]
    new = temp[1] / float(temp[0])
    posit_var[key] = new
                            
#### run psipred to determine secondary structure ####

if exp_seq + '.clean.ss2' not in os.listdir('.'):
    os.system("tcsh ./psipred/BLAST+/runpsipredplus " + exp_seq + ".clean.fasta")
        
## create dictionary of secondary structure elements by position

structure = {}

with open(exp_seq+'.clean.ss2', 'r') as pss:
    l = 0
    for line in pss:
        l = l + 1
        if l > 2 and l < 1002:
            line = re.split(' +', line)
            structure[int(line[1]) - 1] = line[3] # python indexing, starts at 0!
        elif l > 1001:
            line = re.split(' +', line)   
            structure[int(line[0])-1] = line[2]


runcode = ''.join(random.choice(string.ascii_letters) for i in range(6))  # generates a random 6 digit string for labeling
                
##### model section starts here #####

## starts are replications, and begin with the unmutated sequence string

s = 0 # start number
m = 0 # total mutations
rm = 0 # total random mutations

for start in range(starts):
    s = s + 1
    log_out = open(exp_seq + '_' + str(s) + '_' + str(events) + '_' + str(required[4]) + '_' + select + '_' + runcode + '.log','w')
    fasta_out = open(exp_seq + '_' + str(s) + '_' + str(events) + '_' + str(required[4]) + '_' + select + '_' + runcode + '.fasta','w')
    aa_seq_str = aa_seq_str_i
    print('>'+exp_seq+'_'+str(s)+'_'+str(0)+'\n'+aa_seq_str, file = fasta_out)
    aa_seq = list(aa_seq_str)
    n = 0
    temp_out = calculate_pp(aa_seq_str)
    param_output.append([s, n, temp_out[0], temp_out[1], temp_out[2], temp_out[3], temp_out[4]])
    
    ## calculate parameters by 9 aa moving window
    
    gravy = []
    iso = []
    arom = []
    ai = []
    flex = PP.ProteinAnalysis.flexibility(PP.ProteinAnalysis(aa_seq_str))
    
    for i in range(0, len(aa_seq_str)):
        aa_window = aa_seq_str[i:i + 9]
        if len(aa_window) == 9:
            tf, tg, ti, tar, tai = calculate_pp(aa_window)
            gravy.append(tg)
            iso.append(ti)
            arom.append(tar)
            ai.append(tai)
    
    ## calculate mean parameters for each ss element
    
    ss_params = {} ## key = struct, param value = param value}
    
    for struct in ['H', 'E', 'C']:
                   
        tflex = []
        tgravy = []
        tiso = []
        tarom = []
        tai = []
        
        for position in structure.keys():
            if structure[position] == struct:
                if int(position) > 4:
                    try:
                        i = int(position)
                        tflex.append(flex[i])
                        tgravy.append(gravy[i])
                        tiso.append(iso[i])
                        tarom.append(arom[i])
                        tai.append(ai[i])
                    except IndexError:
                        break
                    
        ss_params[struct, 'flex'] = sum(tflex)/len(tflex)
        ss_params[struct, 'gravy'] = sum(tgravy)/len(tgravy)
        ss_params[struct, 'iso'] = sum(tiso)/len(tiso)
        ss_params[struct, 'arom'] = sum(tarom)/len(tarom)
        ss_params[struct, 'ai'] = sum(tai)/len(tai)
            
    coil_param_output.append([s, n, ss_params['C', 'flex'], ss_params['C', 'gravy'], ss_params['C', 'iso'], ss_params['C', 'arom'], ss_params['C', 'ai']])
    beta_param_output.append([s, n, ss_params['E', 'flex'], ss_params['E', 'gravy'], ss_params['E', 'iso'], ss_params['E', 'arom'], ss_params['E', 'ai']])
    helix_param_output.append([s, n, ss_params['H', 'flex'], ss_params['H', 'gravy'], ss_params['H', 'iso'], ss_params['H', 'arom'], ss_params['H', 'ai']])
                
    keep = True ## sets switch to calculate initial probability matrix

    for event in range(events):
        
        n = n + 1
                
        #### Generate mutation ####

        if keep:            
            ## build probability matrix for each position and each possile replacement residue
            probs, sorted_keys, sums = build_mut_matrix(aa_seq_str, aa_freq, sub_matrix, posit_var)

        rp, rr, cr = get_cand(sorted_keys, probs, aa_seq_str, sums)       
                            
        new_aa_seq = aa_seq[:] ## need to copy object, not just assign variable name
        new_aa_seq[rp] = cr                
        new_aa_seq_str = ''.join(new_aa_seq)
                
        #### Evaluate parameters along new sequence ####
        
        ## calculate parameters by 9 aa moving window on the new aa string
    
        gravy = []
        iso = []
        arom = []
        ai = []
        flex = PP.ProteinAnalysis.flexibility(PP.ProteinAnalysis(new_aa_seq_str))
        
        for i in range(0, len(new_aa_seq_str)):
            aa_window = new_aa_seq_str[i:i + 9]
            if len(aa_window) == 9:
                tf, tg, ti, tar, tai = calculate_pp(aa_window)
                gravy.append(tg)
                iso.append(ti)
                arom.append(tar)
                ai.append(tai)
        
        ## calculate mean parameters for each ss element
        
        new_ss_params = {} ## key = struct value = {flex:flex param, ...}
        
        for struct in ['H', 'E', 'C']:
                       
            tflex = []
            tgravy = []
            tiso = []
            tarom = []
            tai = []
            
            for position in structure.keys():
                if structure[position] == struct:
                    if int(position) > 4:
                        try:
                            i = int(position)
                            tflex.append(flex[i])
                            tgravy.append(gravy[i])
                            tiso.append(iso[i])
                            tarom.append(arom[i])
                            tai.append(ai[i])
                        except IndexError:
                            break
                        
            new_ss_params[struct, 'flex'] = sum(tflex)/len(tflex)
            new_ss_params[struct, 'gravy'] = sum(tgravy)/len(tgravy)
            new_ss_params[struct, 'iso'] = sum(tiso)/len(tiso)
            new_ss_params[struct, 'arom'] = sum(tarom)/len(tarom)
            new_ss_params[struct, 'ai'] = sum(tai)/len(tai)
              
        #### this section is adding the selective pressure, if indicated for secondary structure regions

        keep = True
        
        param_str = ''
        
        for direction in director.keys():
            direct_dict = director[direction]
            
            for param in direct_dict.keys():
                old_param = ss_params[direction, param]
                new_param = new_ss_params[direction, param]
                delta_param = float(new_param) - float(old_param)
                param_str = param_str + 'd' + direction + param + ' ' + str(delta_param) + ' '
                
                if direct_dict[param] == 'up':
                    if delta_param <= 0: ## change to < to allow survival of neutral mutations
                        keep = False
                elif direct_dict[param] == 'down':
                    if delta_param >= 0: ## change to > to allow survival of neutral mutations
                        keep = False
                                            
        if keep:
            m = m + 1            
            print(n, m, rm, rp + 1, rr, cr, 'mutation:', param_str)
            print(s, n, m, rm, rp + 1, rr, cr, 'mutation:', param_str , file = log_out)
            pos_w_mut.append(i)
            ss_params = new_ss_params
            aa_seq = new_aa_seq
            aa_seq_str = ''.join(aa_seq)
            
            temp_out = calculate_pp(aa_seq_str)
            param_output.append([s, n, temp_out[0], temp_out[1], temp_out[2], temp_out[3], temp_out[4]])
            
            coil_param_output.append([s, n, ss_params['C', 'flex'], ss_params['C', 'gravy'], ss_params['C', 'iso'], ss_params['C', 'arom'], ss_params['C', 'ai']])
            beta_param_output.append([s, n, ss_params['E', 'flex'], ss_params['E', 'gravy'], ss_params['E', 'iso'], ss_params['E', 'arom'], ss_params['E', 'ai']])
            helix_param_output.append([s, n, ss_params['H', 'flex'], ss_params['H', 'gravy'], ss_params['H', 'iso'], ss_params['H', 'arom'], ss_params['H', 'ai']])

            print('>'+exp_seq+'_'+str(s)+'_'+str(n)+'\n'+aa_seq_str, file = fasta_out)

        else:
            roll = random.sample(dice, 1)
            
            if roll[0] == r_dice:
                m = m + 1
                rm = rm + 1
                print(n, m, rm, rp + 1, rr, cr, 'mutate on dice roll:', param_str)
                print(s, n, m, rm, rp + 1, rr, cr, 'mutate on dice roll:', param_str ,file = log_out)
                pos_w_mut.append(i)
                aa_seq = new_aa_seq
                aa_seq_str = ''.join(aa_seq)
                
                temp_out = calculate_pp(aa_seq_str)
                param_output.append([s, n, temp_out[0], temp_out[1], temp_out[2], temp_out[3], temp_out[4]])
            
                coil_param_output.append([s, n, ss_params['C', 'flex'], ss_params['C', 'gravy'], ss_params['C', 'iso'], ss_params['C', 'arom'], ss_params['C', 'ai']])
                beta_param_output.append([s, n, ss_params['E', 'flex'], ss_params['E', 'gravy'], ss_params['E', 'iso'], ss_params['E', 'arom'], ss_params['E', 'ai']])
                helix_param_output.append([s, n, ss_params['H', 'flex'], ss_params['H', 'gravy'], ss_params['H', 'iso'], ss_params['H', 'arom'], ss_params['H', 'ai']])
                
                print('>'+exp_seq+'_'+str(s)+'_'+str(n)+'\n'+aa_seq_str, file = fasta_out)
            else:
                print(n, m, rm, rp + 1, rr, cr, 'none:', param_str)
                print(s, n, m, rm, rp + 1, rr, cr, 'none:', param_str, file = log_out)
               
    fasta_out.close()
    
    #### model over at this point, rest is just printing output ####
    
    ## tally the number of mutations for each position
    ## positions numbered according to python indexing
    ## minimum number of mutations is 1 - meaning only one residue has occupied that position    


    with open(exp_seq + '_' + str(s) + '_' + str(events) + '_' + str(required[4]) + '_' + select + '_' + runcode + '.fasta', 'r') as fasta_hist:
        position_dict = {} ## dictionary holding tally of difference bases for each position
        for line in fasta_hist:
            line = line.rstrip('\n')
            if line.startswith('>') == False:
                for r,residue in enumerate(line):
                    if r not in position_dict.keys():
                        position_dict[r] = residue
                    elif residue != position_dict[r][-1:]:
                        position_dict[r] = position_dict[r]+residue
                        
    ## print the number of mutations and the sequence of mutations for each position
    
    mutation_tally = open(exp_seq + '_' + str(s) + '_' + str(events) + '_' + str(required[4]) + '_' + select + '_' + '_tally' + '_' + runcode + '.txt','w')
    for key in position_dict.keys():
        tally = len(position_dict[key])
        print(str(key)+'\t'+structure[key]+'\t'+str(tally)+'\t'+position_dict[key] , file = mutation_tally)
    mutation_tally.close()
    
print('mutations =', m, 'random mutations =', rm, 'events =', events * starts)


with open(exp_seq + '_' + str(starts) + '_' + str(events) + '_' + str(required[4]) + '_' + select + '_' + runcode + '.txt','w') as param_file:
    print('run' + '\t' + 'generation' + '\t' + 'flex' + '\t' + 'gravy' + '\t' + 'iso' + '\t' + 'arom' + '\t' + 'ai', file = param_file)    
    nested_csv(param_file, param_output)
    
with open(exp_seq + '_coil_' + '_' +  str(starts) + '_' + str(events) + '_' + str(required[4]) + '_' + select + '_' + runcode + '.txt','w') as param_file:
    print('run' + '\t' + 'generation' + '\t' + 'flex' + '\t' + 'gravy' + '\t' + 'iso' + '\t' + 'arom' + '\t' + 'ai', file = param_file)    
    nested_csv(param_file, coil_param_output)
    
with open(exp_seq + '_beta_' + '_' + str(starts) + '_' + str(events) + '_' + str(required[4]) + '_' + select + '_' + runcode + '.txt','w') as param_file:
    print('run' + '\t' + 'generation' + '\t' + 'flex' + '\t' + 'gravy' + '\t' + 'iso' + '\t' + 'arom' + '\t' + 'ai', file = param_file)    
    nested_csv(param_file, beta_param_output)
    
with open(exp_seq + '_helix_' + '_' + str(starts) + '_' + str(events) + '_' + str(required[4]) + '_' + select + '_' + runcode + '.txt','w') as param_file:
    print('run' + '\t' + 'generation' + '\t' + 'flex' + '\t' + 'gravy' + '\t' + 'iso' + '\t' + 'arom' + '\t' + 'ai', file = param_file)    
    nested_csv(param_file, helix_param_output)
    
with open(exp_seq + '_posit_var_' + runcode + '.txt', 'w') as var_out:
    for key in posit_var.keys():
        print(key, posit_var[key], file = var_out) 
    
print('Run ID: ' + runcode)

log_out.close()
