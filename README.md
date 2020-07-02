# O157 Lineage Assignment
A script that returns the LSPA6 lineage result for O157 Escherichia coli when given genomic data in a FASTA format.

The LSPA6 lineage typing is based on: 




Assumptions
If you are a USDA user with access, this program runs natively on CERES as of 7/1/2020.
Python is installed on machine. Download python 2.7.
User has bash terminal.

Input
Open Terminal
Call program ./LSPA6Long.sh [output CSV File Name] [Fasta File(s)]

Output
A .csv file will be created containing concise information about each strain or contig (as separated by >).
Unless directed to a different location, the output CSV file will be created wherever the script is located.

E.g.: Output
Verbose output on the screen is also an output. This output can be sent to a file. ./LSPA6Long.sh [output CSV File Name] [Fasta File(s)] > Information.txt
