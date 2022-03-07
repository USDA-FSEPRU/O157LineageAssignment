#!/usr/bin/env python2.7

"""
This script assigns LSPA6 status to an inputed FASTA file.
The scheme is used to infer genetic distance between E. coli strains. 
More information about the schemes are available: 


"""
import re
import sys
import csv
import math

# Check to see if a value is a float
def isfloat(value):
  try:
    float(value)
    return True
  except:
    return False

inFiles = sys.argv[2:]
csvFile= sys.argv[1]
start=2 

count = len(sys.argv)



# Primers to test
#folD-sfmA
folD1 = 'TACGTAGGTCGAAGGG' 
folD1bR = 'CCAGATTTACAACGCC'  
folD1R = 'CCCTTCGACCTACGTA'
folD1b = 'GGCGTTGTAAATCTGG'

#X5935
z59351 = 'GTGTTCCCGGTATTTG'
z59351bR = 'CTCACTGGCGTAACCT'
z59351R = 'CAAATACCGGGAACAC'
z59351b = 'AGGTTACGCCAGTGAG'

# yhcG
yhcG1 = 'CTCTGCAAAAAACTTACGCC'
yhcG1bR = 'CAGGTGGTTGATCAGCG'
yhcG1R= 'GGCGTAAGTTTTTTGCAGAG'
yhcG1b= 'CGCTGATCAACCACCTG'

#rbsB1
rbsB1 = 'AGTTTAATGTTCTTGCCAGCC'
rbsB1bR = 'ATTCACCGCTTTTTCGCC'
rbsB1R = 'GGCTGGCAAGAACATTAAACT'
rbsB1b = 'GGCGAAAAAGCGGTGAAT'

        
#rtcB
rtcB1 = 'GCGCCAGATCGATAAAGTAAG'
rtcB1bR ='GCCGTTGTAAACGTGATAAAG'
rtcB1R = 'CTTACTTTATCGATCTGGCGC'
rtcB1b = 'CTTTATCACGTTTACAACGGC'
        
#arp-iclR
arp1 = 'GCTCAATCTCATAATGCAGCC'
arp1bR = 'CACGTATTACCGATGACCG'
arp1R = 'GGCTGCATTATGAGATTGAGC'
arp1b= 'CGGTCATCGGTAATACGTG' 

sequence = ''
genome= ''

def near_hit(primer, str1):
    if len(primer) != len(str1): #rare pathological case of partial match near ends of sequences
        return False 
    miss=0 #miss has not occurred
    for i in range(0,len(primer)):
        if str1[i] != primer[i]:
            if miss>1: #miss already happened
                return False #two misses
            miss=miss+1 #miss has occurred
    return True #zero or one misses

def frag_list_check(index_list1, index_list2, min_frag_length, max_frag_length):
    assert min_frag_length < max_frag_length
    if (not index_list1) or (not index_list2):
        return False
    for item1 in index_list1:
        for item2 in index_list2:
            if max_frag_length>abs(item1-item2)>min_frag_length:
                return True
    return False

def primer_search(sequence, primer):
    l=len(primer)
    #break primer into two halves and make empty lists to store matches: 
    front=primer[:l/2]
    front_list=[]
    front_list_ind=[]
    back=primer[l/2:]
    back_list=[]
    back_list_ind=[]

    #search genome twice. This is far more efficient than searching many times.
    ind=-1;
    while ind < len(sequence):
        ind=sequence.find(front,ind+1)
        if ind <= -1: 
            break
        front_list_ind.append(ind);
        front_list.append(sequence[ind:ind+len(primer)])
    ind=-1;
    while ind < len(sequence):
        ind=sequence.find(back,ind+1)
        if ind <= -1: 
            break
        back_list_ind.append(ind-l/2);
        back_list.append(sequence[ind-l/2:ind+len(primer)-l/2])    
    #merge lists
    full_list=front_list+back_list
    full_list_ind=front_list_ind+back_list_ind

    #first look for exact hits.
    hit_list_ind=[]
    i=0
    for item in full_list:
        if item == primer:
            hit_list_ind.append(full_list_ind[i])
        i=i+1
    if hit_list_ind: #if non-empty 
        hit_list_ind=list(set(hit_list_ind))
        print "Primer: ", primer, 'Exact Match: ', hit_list_ind
        return hit_list_ind #return list of full primer match indices

    #now look for "near hits". The near_hit function determines what a near hit is. 
    hit_list_ind=[]
    i=0;
    for item in full_list:
        if near_hit(primer, item):
            hit_list_ind.append(full_list_ind[i])
        i=i+1

    if hit_list_ind:
        print "Primer: ", primer, "Near Match: ", hit_list_ind
    else:
        print "Primer: ", primer, "No primer match."
    return hit_list_ind




#Create a .csv file
with open(csvFile, 'wb') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(['A verbose .txt file containing gene fragments and location can also be generated. See documentation for more information.'])
    filewriter.writerow(['Accession', 'foldD-sfmA', 'Z5935', 'yhcG', 'rtcB', 'rbsB', 'arp-iclR', 'LSPA type', 'lineage'])
#Loop through files and process them
    for i in range(start, count):
        inFile=sys.argv[i]
        i = i + 1
        accession = sys.argv[i]
        # this next line is a problem:
        # Assumes the fasta headers have a description.  If this is not true then the 'accession' variable becomes the whole sequence
        with open(inFile, 'r') as myfile:
            inFile=myfile.read().replace('\n', '')
        
        sequence = inFile
        #accession = sequence[1:].partition(' ')[0] #was 1
        sequence=sequence.upper()
        print 'Accession: '+accession
        
        
        print '\n----------Gene Results-----------'
        
        # SEARCH for folD
        print '\nfolD-sfmA results:'
        folD = 0
        
        temp1 = primer_search(sequence, folD1)
        temp2 = primer_search(sequence, folD1b)
        temp3 = primer_search(sequence, folD1R) 
        temp4 = primer_search(sequence, folD1bR)

        if frag_list_check(temp1, temp2, 144, 146):
            folD = 1
            print 'foldD=1'
        elif frag_list_check(temp3, temp4, 144, 146):
            folD = 1
            print 'foldD=1'
        elif frag_list_check(temp1, temp2, 153, 155):
        	folD = 2
        	print 'foldD=2'
        elif frag_list_check(temp3, temp4, 153, 155):
        	folD = 2
        	print 'foldD=2'
        else:
            folD = 0
            print 'foldD=0'
        
        
        
        # SEARCH for z5935
        print '\nZ5935 results:'
        Z5935 = 0

        temp1 = primer_search(sequence, z59351)
        temp2 = primer_search(sequence, z59351b)
        temp3 = primer_search(sequence, z59351R) 
        temp4 = primer_search(sequence, z59351bR)
        
        if  frag_list_check(temp1, temp2, 116, 118):
            Z5935 = 1
            print 'Z5935=1'
        elif frag_list_check(temp3, temp4, 116, 118):
            Z5935 = 1
            print 'Z5935=1'
        elif  frag_list_check(temp1, temp2, 125, 127):
            Z5935 = 2
            print 'Z5935=2'
        elif frag_list_check(temp3, temp4, 125, 127):
            Z5935 = 2
            print 'Z5935=2'
        #other fragments that would appear in gel
        elif  frag_list_check(temp1, temp2, 119, 124):
            Z5935 = 3
            print 'Z5935=3'
        elif frag_list_check(temp3, temp4, 119, 124):
            Z5935 = 3
            print 'Z5935=3'
        elif  frag_list_check(temp1, temp2, 128, 500):
            Z5935 = 3
            print 'Z5935=3'
        elif frag_list_check(temp3, temp4, 128, 500):
            Z5935 = 3
            print 'Z5935=3'
        else:
            Z5935 = 0
            print 'Z5935=0'

        
        # SEARCH for yhcG
        print '\nyhcG results:'
        yhcG = 0

        temp1 = primer_search(sequence, yhcG1)
        temp2 = primer_search(sequence, yhcG1b)
        temp3 = primer_search(sequence, yhcG1R) 
        temp4 = primer_search(sequence, yhcG1bR)
        
        #allele 1
        if frag_list_check(temp1, temp2, 373, 378):
            yhcG = 1
            print 'yhcG=1'
        elif frag_list_check(temp3, temp4, 373, 378):
            yhcG = 1
            print 'yhcG=1'
        #allele 2
        elif frag_list_check(temp1, temp2, 451, 456):
            yhcG = 2
            print 'yhcG=2'
        elif frag_list_check(temp3, temp4, 451, 456):
            yhcG = 2
            print 'yhcG=2'
        # other alleles
        elif frag_list_check(temp1, temp2, 194, 372):
            yhcG = 3
            print 'yhcG=3'
        elif frag_list_check(temp3, temp4, 194, 372):
            yhcG = 3
            print 'yhcG=3'
        elif frag_list_check(temp1, temp2, 379, 450):
            yhcG = 3
            print 'yhcG=3'
        elif frag_list_check(temp3, temp4, 379, 450):
            yhcG = 3
            print 'yhcG=3'
        else:
            yhcG = 0
            print 'yhcG=0'
            
        # SEARCH for rbsB
        print '\nrbsB results:'
        rbsB = 0

        temp1 = primer_search(sequence, rbsB1)
        temp2 = primer_search(sequence, rbsB1b)
        temp3 = primer_search(sequence, rbsB1R) 
        temp4 = primer_search(sequence, rbsB1bR)
        
        #allele 1
        if frag_list_check(temp1, temp2, 196, 201):
            rbsB = 1
            print 'rbsB=1'
        elif frag_list_check(temp3, temp4, 196, 201):
            rbsB = 1
            print 'rbsB=1'
        #allele 2
        elif frag_list_check(temp1, temp2, 187, 192):
            rbsB = 2
            print 'rbsB=2'
        elif frag_list_check(temp3, temp4, 187, 192):
            rbsB = 2
            print 'rbsB=2'
        # other alleles
        elif frag_list_check(temp1, temp2, 193, 195):
            rbsB = 3
            print 'rbsB=3'
        elif frag_list_check(temp3, temp4, 193, 195):
            rbsB = 3
            print 'rbsB=3'
        elif frag_list_check(temp1, temp2, 202, 500):
            rbsB = 3
            print 'rbsB=3'
        elif frag_list_check(temp3, temp4, 202, 500):
            rbsB = 3
            print 'rbsB=3'
        else:
            print 'rbsB=0'
        
        
        # SEARCH for rtcB
        print '\nrtcB results:'
        rtcB = 0

        temp1 = primer_search(sequence, rtcB1)
        temp2 = primer_search(sequence, rtcB1b)
        temp3 = primer_search(sequence, rtcB1R) 
        temp4 = primer_search(sequence, rtcB1bR)
        
        #allele 1
        if frag_list_check(temp1, temp2, 248, 250):
            rtcB = 1
            print 'rtcB=1'
        elif frag_list_check(temp3, temp4, 248, 250):
            rtcB = 1
            print 'rtcB=1'
        #allele 2
        elif frag_list_check(temp1, temp2, 257, 259):
            rtcB = 2
            print 'rtcB=2'
        elif frag_list_check(temp3, temp4, 257, 259):
            rtcB = 2
            print 'rtcB=2'
        # other alleles
        elif frag_list_check(temp1, temp2, 251, 256):
            rtcB = 3
            print 'rtcB=3'
        elif frag_list_check(temp3, temp4, 251, 256):
            rtcB = 3
            print 'rtcB=3'
        elif frag_list_check(temp1, temp2, 260, 500):
            rtcB = 3
            print 'rtcB=3'
        elif frag_list_check(temp3, temp4, 260, 500):
            rtcB = 3
            print 'rtcB=3'
    	elif frag_list_check(temp1, temp2, 200, 247):
            rtcB = 3
            print 'rtcB=3'
        elif frag_list_check(temp3, temp4, 200, 247):
            rtcB = 3
            print 'rtcB=3'
        else:
            rtcB = 0
            print 'rtcB=0'


        # SEARCH for arp-iclR
        print '\narp-iclR results:'
        arpiclR = 0

        temp1 = primer_search(sequence, arp1)
        temp2 = primer_search(sequence, arp1b)
        temp3 = primer_search(sequence, arp1R) 
        temp4 = primer_search(sequence, arp1bR)
        
        #allele 1
        if frag_list_check(temp1, temp2, 292, 296):
            arpiclR = 1
            print 'arp=1'
        elif frag_list_check(temp3, temp4, 292, 296):
            arpiclR = 1
            print 'arp=1'
        #allele 2
        elif frag_list_check(temp1, temp2, 310, 314): #was 312
            arpiclR = 2
            print 'arp=2'
        elif frag_list_check(temp3, temp4, 310, 314): #was 312
            arpiclR = 2
            print 'arp=2'
        # allele 3
        elif frag_list_check(temp1, temp2, 297, 310): #was 311
            arpiclR = 3
            print 'arp=3'
        elif frag_list_check(temp3, temp4, 297, 310): #was 311
            arpiclR = 3
            print 'arp=3'
        elif frag_list_check(temp1, temp2, 200, 291):
            arpiclR = 3
            print 'arp=3'
        elif frag_list_check(temp3, temp4, 200, 291):
            arpiclR = 3
            print 'arp=3'
    	elif frag_list_check(temp1, temp2, 315, 500):
            arpiclR = 3
            print 'arp=3'
        elif frag_list_check(temp3, temp4, 315, 500):
            arpiclR = 3
            print 'arp=3'
        else:
            print 'arp=0'
        

        oldGroup=''
        newGroup=''
        warning=''
        print '\n-----------CONCLUSIONS----------'

    	lineage='check genotype'
        	
        concat = str(folD) + str(Z5935) + str(yhcG) + str(rtcB) + str(rbsB) + str(arpiclR)
        print 'LSPA type: ' + concat
        
        if concat == '111111':
    		lineage = 'LI'
    		print "lineage: " + lineage
    	elif concat == '211111':
    		lineage = 'LI/II'
    		print "lineage: " + lineage
    	elif concat == '222222' or concat == '222212' or concat == '222221' or concat == '222211':
    		lineage = 'LII'
    		print "lineage: " + lineage
    	elif concat == '000000':
    		lineage = "check sequencing data"
    	else:
    		lineage= 'manually assign'
    		print "lineage: " + lineage
        
        filewriter.writerow([accession,folD,Z5935,yhcG,rbsB,rtcB,arpiclR,concat,lineage])
        
    
        print '\nHave a great day!\n\n'   

csvfile.close()


       
