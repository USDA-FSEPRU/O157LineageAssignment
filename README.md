# O157 Lineage Assignment
A script that returns the Lineage-Specific Polymorphism Assay-6 (LSPA-6) result for O157 *Escherichia coli* when given genomic data in a FASTA format.

The LSPA6 lineage typing is based on: 

[Yang Z, Kovar J, Kim J, Nietfeldt J, Smith DR, Moxley RA, Olson ME, Fey PD, Benson AK. 2004. Identification of common subpopulations of non-sorbitol-fermenting, beta-glucuronidase-negative *Escherichia coli* O157:H7 from bovine production environments and human clinical samples. Appl Environ Microbiol 70:6846-54.](https://aem.asm.org/content/70/11/6846/article-info)

<br>

### Assumptions
* If you are a USDA user with access, this program runs natively on CERES as of 7/1/2020.
* Python is installed on machine.
* User has bash terminal.

<br>

### Input
* Open Terminal
* Call program 

```./LSPA6Long.sh [output CSV File Name] [Fasta File(s)]```

<br>

### Output
* A .csv file will be created containing concise information about each strain or contig (as separated by >).
* Unless directed to a different location, the output CSV file will be created wherever the script is located.

![Standard Out](https://github.com/nielsend/O157LineageAssignment/blob/master/CSVImage.png)


* Verbose output on the screen is also an output. This output can be sent to a file. 

```./LSPA6Long.sh [output CSV File Name] [Fasta File(s)] > Information.txt```

![Standard Out](https://github.com/nielsend/O157LineageAssignment/blob/master/StandardOut.png)


* Please note that LSPA-6 results other than 111111, 211111, and 222222 will be assigned the lineage of "manually assign". Please examine these genomes with greater scrutiny.
