## methylsanger-matlab

Make the command file to feed to MATLAB:

To make running the bisulfitePCRSeqAnalysis.m script easier, we have created the helper python script, writeMfile.py. This python script generates the MATLAB command file for you. 
To run this script, you will need to give the name of the command file, the name of the scf file, the name of the alignment file (generated using bsAlign.py and bsDraw.py) and the name of the output image file (use .png extension). After running it should print the coordinates of the beginning and ending bounds that will be analyzed. An example run is below:

```
writeMfile.py cmd.m exp-GRP78RC.scf a-nd-GRP78_ref.fa exp-GRP78RC.png
```

Example output to screen will be the beginning and end of the sequnce to be analyzed, e.g., Bounds found at positions 80, 442.

This produces the file cmd.m with the following format (NB: This should all be on one line!!!; I've formatted to fit in display!!!):

```
bisulfitePCRSeqAnalysis('exp-GRP78RC.scf', 'AAACACCCCAATAGGTCAATCTGTCTGTGCTGTCTTGGCCGGCGTCGACCTCACCGTCGCCTACTCGGCTTATATACCCTCCCCCAGCCCCGTCGTGGAGGCCGCC
GATTGGTGAAGGCCCTGCTCGTTGGAGGCCGTTCATTGGCCCAGGCCACCAAGCTGGCCGCCGCCGATTCGAAGCGGCCCCCTCCGCAATAAACGTCACTGCTCC
GCCCCCCCCAGCTGATTCATTGGCTGCTATTCGTTTCTAACGTTCACCAATGGTAGATAACATCCGCCCCATCCGCTGGTCCCATTGGTTCAGC
AGCTGTCTATCTCTCCTGCGACTTCTGACCCCGAGGCATTTCCGCTGGTAACCGCACAAGTCCCGCCTTCACTCCCGGCCCTAGGGGGTCGG
AGTAGGTCCAGCAGGAGTGACCCCCCGGGGCTGTGGGATCTGAAACTTTTTCTTCTC', 'CCCCCAGCC', 'CTTTTTCTT', 'NCCCCCAGN', 'CTTTTTCTN', 'exp-GRP78RC.png');
```
Run the MATLAB  Script:

matlab cmd.m
