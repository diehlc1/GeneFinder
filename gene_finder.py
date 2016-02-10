# -*- coding: utf-8 -*-
"""
DNA sequencing and amino acid identification
Takes in a string of DNA

@author: Cecilia Diehl

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('F') #adding this to make sure the function finds non nucleotide letters
    'Not a base pair'
    >>> get_complement('G')
    'C'
    """
    #replacing nucleotides with their correct base pair
    if nucleotide == 'A':
    	complement = 'T'
    elif nucleotide == 'C':
        complement = 'G'
    elif nucleotide == 'T':
        complement = 'A'
    elif nucleotide == 'G':
        complement = 'C'
    else :
    	return "Not a base pair"

    return complement


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("A") #does the function get caught up on short strings? nope, cool
    'T'
    """
    full_complement = "" #making an empty string
    reverse_complement = ""
    i = 0
    for letter in dna:
        complement = get_complement(letter)
        full_complement = full_complement + complement #filling full_complement with results from get_complement
    
    for i in range(1, len(full_complement)+1):
        nucleotide = full_complement[-i]
        reverse_complement = reverse_complement + nucleotide
    
    return reverse_complement #if I gave it 'ACT' it gives me 'TGA'


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA") #testing TGA stop codon
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG") #testing TAG stop codon
    'ATGAGA'
    >>> rest_of_ORF("ATGGCATCTTAA") #testing TAA stop codon
    'ATGGCATCT'
    >>> rest_of_ORF("ATGATCTTAA") #tsting if the stop codon is not in a set of three, returns dna
    'ATGATCTTAA'
    """
    stop = ["TAA", "TAG", "TGA"]
    for i in range(0, len(dna),3):  
        if dna[i: i+3] in stop:#determining if the select 3 values are a stop codon
            return dna[:i]
    return dna 


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGGGATGGTAG") #if there is a nested orf it is not returned individualy
    ['ATGGGATGG']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCTAA") #if the stop codon is not in a set of three
    ['ATGCATGAATGTAGA', 'ATGTGCCTAA']
    >>> find_all_ORFs_oneframe("TGGATGCATGAATGTAGATAGATGTGCCC") #if the start codon is not at the begining
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    start = ["ATG"]
    i = 0
    orf_list = [] #creating an empty list to fill later
    while i < len(dna):
    #for i in range(0, len(dna)):
        if dna[i:i+3] in start:
            orf_list.append(rest_of_ORF(dna[i:]))
            i += len(rest_of_ORF(dna[i:]))
        i+=3
    return orf_list


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs_oneframe("TATATGCATGAATGGTGAA") #if the start codon is offset 
    ['ATGCATGAATGG']
    >>> find_all_ORFs_oneframe("ATGCATCGCTGAATGAATTGACC") #more than one orf in a frame, nested orf ignored
    ['ATGCATCGC', 'ATGAAT']
    """
    orf_list = []
    orf_list.extend(find_all_ORFs_oneframe(dna[:]))  #starting find_all_ORFs_oneframe offset each time 
    orf_list.extend(find_all_ORFs_oneframe(dna[1:])) #extend replaces append for these
    orf_list.extend(find_all_ORFs_oneframe(dna[2:])) #when using append it makes a list of a list  
    return orf_list


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("ATGCATGAATGTAG") #dna from find_all_ORFs() should have the same orfs plus more from the reverse_complement dna 
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG', 'ATGCAT']
    """
    orf_list = []
    orf_list.extend(find_all_ORFs(dna)) #putting some orfs in a list
    orf_list.extend(find_all_ORFs(get_reverse_complement(dna))) #including the compliment of the dna strand
    return orf_list[::]

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string. I think that this test is already pretty good because it is working
        from functions that have already been thoughrly tested, so as long as it can pick out
        the longest orf string in the given list it works
    >>> longest_ORF("ATGCGAATGTAGCATCAAA") #longest orf is last in the list
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGCATGAATGTAG") #longest orf is first in the list
    'ATGCATGAATGTAG'
    """
    orf_list = []
    orf_list.extend(find_all_ORFs_both_strands(dna)) #including all orfs
    return max(orf_list, key=len) #taking the orf with the max length


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    length = 0
    for trial in range(num_trials):
        string = shuffle_string(dna)
        if len(str(longest_ORF(string)))>length: 
            length = len(str(longest_ORF(string))) #replaces the value (orf) if it finds one longer
    return length
    


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
                 #the start codon is offset in the seccond test, however this isn't a problem because the 
                 #dna used in the final run will have already been separated into orfs starting with start codons

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT") 
        'MPA'
        >>> coding_strand_to_AA("GGCATGCCCGCTTT") #the start codon is offset 
        'GMPA'
        """
    sequence = ""
    for i in range(0,len(dna),3):
        codon = dna[i:i+3]
        if len(codon) == 3: 
            amino_acid = aa_table[codon]
            sequence = sequence + amino_acid
    return sequence 




def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    sequence_list = []
    threshold = longest_ORF_noncoding(dna, 1500) #returns a length, int
    orfs = find_all_ORFs_both_strands(dna) #returns list
    for orf in orfs:
        if len(orf) >= threshold:
            amino_acid = coding_strand_to_AA(orf)
            sequence_list.append(amino_acid)
    return sequence_list
    

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    #dna = load_seq("./data/X73525.fa") #uncomment to run w/ imported dna
    #print gene_finder(dna) #uncomment to print

