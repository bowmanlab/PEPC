# PEPC
The Protein Evolution Parameter Calculator

Created on Wed Dec 05 15:04:04 2012\
Jeff Bowman, bowmanjs@gmail.com

## Recommended Settings
This guide assumes you are either running a Windows or Linux distribution. No guide is yet avaliable for MacOS. 

If you are running windows, I recommend using the lastest version of [Ubuntu](https://ubuntu.com/) for WSL. In addition, the psipred package, which you will download later, requires use of the Tenex C Shell (AKA `tcsh`). You can follow [this guide](https://randomknowhow.com/tech/infrastructure/devops/tutorials/c-shell-tc-shell-install-for-linux-windows-mac/) to install and configure both Ubuntu and tcsh. Do not follow the directions to make bash automatically launch tcsh: stop after step 4. You might as well download the "C shell" (AKA csh), the same way. **runpsipredplus** seems to work with both. 
 
PEPC is currently written for Python 2.7 and Python 3.7+. It is recommend you use Python 3 to run PEPC. 

Make sure to free up a few GB of storage for the Pfam-A protein database you will need to download. 

## Python Packages
You will need to install [biopython](https://github.com/biopython/biopython). Note that Python 2 references the biopython package as `Bio`, so you may need to rename the folder **biopython** to **Bio** in wherever you store your Python 2 library site-packages. I reccommend configuring and testing according to the biopython manual. PEPC uses the **Bio.SubsMat.MatrixInfo** and **Bio.SeqaUtils.ProtParam** methods. 

## BLAST dependencies

In addition to Biopython PEPC requires the use of psipred.  You will need to
download and install this dependency separately, along with your preferred
protein database for blasting (recommend Pfam-A.fasta). You can clone the psipred files from GitHub [here](https://github.com/psipred/psipred). In bash, simply enter the command `tcsh` to enter the "tee-shell", and follow the direction outlined in the psipred README. The **runpsipredplus** script uses NCBI BLAST+, which you can install from [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). We will need the **psiblast.exe** executable from this download. 

Specifically PEPC makes use of the **runpsipredplus** script, so make sure that this script is properly configured and available in your path. To do so, edit the **runpsipredplus** as follows:

Example
```
# Where the NCBI BLAST+ programs have been installed
set ncbidir = ./psipred/BLAST+
 
# Where the PSIPRED V3 programs have been installed
set execdir = ./psipred/bin
 
# Where the PSIPRED V3 data files have been installed
set datadir = ./psipred/data
```

The period in `./filename` represents the current directory, which will be the PEPC parent directory when running the model. Make sure to set `ncbidir` as the directory where your **psiblast.exe** executable lives. 

I recommend running the PSI-Blast against Pfam-A.fasta, after conversion this file to a blast database. Go to <ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release> and click on 'Pfam-A.fasta.gz' to download. Alternatively, you can use the command `wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release`, but downloading this way was taking too long for me. This site has many other databases, such as unitprot, if you plan on using a different protein database. Move this file to wherever you want. I just put in in BLAST+ in the psipred folder, since that is what where the **runpsipredplus** script is, then cd to that directory. Unzip the file using 
`gunzip -v Pfam-A.fasta.gz`. Then, you can make this into a BLAST database using the command `makeblastdp -dbtype prot -in Pfam-A.fasta -out Pfam-A.fasta`. This will create a bunch of .phr, .pin, and .psq files, and one .pal file. Go ahead and move all of these to the parent PEPC directory. This looks messy and can probably be fixed, but for some reason this is the only way I could make the **runpsipredplus** script play nice with the protein database we just made. 

You will also need to add **runpsipredplus** as as a command to your source file (.bashrc or .bash_profile). In you .bashrc file write
```
alias runpsipredplus='PATH_TO_COMMAND'
```
This is to help it work with the tcsh as we are jumping directories. 

Next, modify the blast line of the runpsipredplus file like so:
```
# The name of the BLAST+ data bank
set dbname = Pfam-A.fasta
```
or whatever database you are using. 

Now, your psipred codebase should be properly configured to work with the PEPC model. At this point I recommend running `python PEPC_vX.py` to check for errors. If at some point you see the error `FATAL: Error whilst running blastpgp - script terminated!`, this should actually say `psiblast` instead of `blastpgp`. You may want to change this yourself to avoid a future headache. 

You may also wish to specify the number of cpus for blast to use with the num_threads parameter on this line.  Refer to the blast documentation at
<http://www.ncbi.nlm.nih.gov/books/NBK1763> for more details.


To determine site-specific mutation rates PEPC relies on blastp.  It is
recommended that you run blastp against the Pfam-A database used by psipred, but you need to tell PEPC where to look for this database by setting the
`db` parameter at the top of the PEPC python file, like so:

```
#### user setable parameters ####
 
db = '/psipred/BLAST+/Pfam-A' ## database for blastp search
```

After this point, you should have all the depenecies you need to run PEPC on your protein sequence. 

## Usage
Run PEPC as\
`python PEPC_vX.py [options] [fasta] [starts] [events] [dice]`

where:\
`X` = the version number\
`[options]` = a combination of secondary structure elements coil (C), beta sheets
(E), or alpha helices (H) and parameters (flex, gravy, iso, arom, ai) and
direction (up, down), e.g. `C:flex=up,gravy=down`\
`[fasta]` = a protein fasta in the working directory containing a single sequence, omit the .fasta extension.  The file name can NOT have a colon in it. Sorry. NEW: your `[fasta]` can be a path to a fasta file (still with no .fasta extension). Please try to avoid the use of spaces in your higher directories as PSIPRED doesn't seem to handle this well. 
`[starts]` = the number of replicate runs\
`[events]` = the number of mutations you would like to attempt to see an improvement in the parameter space.  Setting to 0 will calculate the protein parameters by secondary structure.\
`[dice]` = the number of sides to the dice rolled to determine if a mutation that does not improve fitness should survive.  0 turns off dice, meaning only mutations that improve fitness will be kept.

PEPC outputs four sets of files.  For each start a new fasta is created containing the mutated aa sequence. For each program run an output file containing the protein parameters for each mutation is created. This output file is in the order: 

start, event (0 = the unmutated sequence), isoelectric point, flexibility, aromaticity, gravy, aliphatic index

The tally.txt file contains the mutations by position. A log file is also generated for each program run. The line is printed for each event in the order start, event, number of mutations, number of random mutations, position of possible mutation, current residue at that position, candidate residue, success of mutation, and resulting change to the indicated parameters.

Protein parameter isoelectric point, flexibility, aromaticity, and gravy are calculated using the ProtParam.ProteinAnalysis module in Biopython.
Aliphatic index is calculated using the following equation:\

100*[A] + 2.9*100*[V] + 3.9*100*[I]+100*[L] according to Akai, 1980.


Full command example:\
`python PEPC_v0.17.py C:flex=down,gravy=up E:arom=up fasta_file_basename 3 1000 10000`



## Version updates 

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
reassignment of the aa_seq object (see <http://stackoverflow.com/questions/19341365/setting-two-arrays-equal> ).

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
