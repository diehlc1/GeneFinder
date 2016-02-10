# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

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
    """
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

#dna = "ACTG"

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    full_complement = "" #saying that reverse_complement is a string
    reverse_complement = ""
    i = 0
    for letter in dna:
        complement = get_complement(letter)
        full_complement = full_complement + complement #filling reverse_complement with get_complement   print reverse_complement
    
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
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGGAGATAGG")
    'ATGGAGA'
    """
    stop = ["TAA", "TAG", "TGA"]
    start = ["ATG"]
    i = 0
    orf = ""
    orf_list = []
    for i in range(0, len(dna),3):
        if dna[i: i+3] is in stop:
            orf_list.append(dna[:i])
    # s1 = "TAA"
    # s2 = "TAG"
    # s3 = "TGA"
    # i = 0
    # orf = ""
    # for i in range(0, len(dna)):
    #     orf = orf + dna[i]
    #     if (dna[i-2]+dna[i-1]+dna[i]) == s1:
    #         orf = orf[:-3]
    #         return orf
    #     elif (dna[i-2]+dna[i-1]+dna[i]) == s2:
    #         orf = orf[:-3]
    #         return orf
    #     elif (dna[i-2]+dna[i-1]+dna[i]) == s3:
    #         orf = orf[:-3]
    #         return orf


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
    """

    s1 = "TAA"
    s2 = "TAG"
    s3 = "TGA"
    stop = ["TAA", "TAG", "TGA"]
    start = ["ATG"]
    i = 0
    orf = ""
    orf_list = []
    for i in range(0, len(dna),3):
        if dna[i: i+3] is in stop:
            orf_list.append(dna[:i])

        # orf = orf + dna[i]
        # if len(orf_list) == 0:
        #     if (dna[i-2]+dna[i-1]+dna[i]) == s1:
        #         orf_list.append(orf[:-3])
        #         orf = ""
        #     elif (dna[i-2]+dna[i-1]+dna[i]) == s2:
        #         orf_list.append(orf[:-3])
        #         orf = ""
        #     elif (dna[i-2]+dna[i-1]+dna[i]) == s3:
        #         orf_list.append(orf[:-3])
        #         orf = ""
        # elif len(orf_list) != 0:
        #     if (dna[i-2]+dna[i-1]+dna[i]) == s1:
        #         orf_list.append(orf[:-3])
        #         orf = ""
        #     elif (dna[i-2]+dna[i-1]+dna[i]) == s2:
        #         orf_list.append(orf[:-3])    
        #         orf = ""
        #     elif (dna[i-2]+dna[i-1]+dna[i]) == s3:
        #         orf_list.append(orf[:-3])
        #         orf = ""
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
    """
    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
