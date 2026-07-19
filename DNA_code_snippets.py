"""

    code pour solutionner les problèmes proposés par Rosalind.info/problems
    
    https://rosalind.info/problems

"""

DNA_NUCLEOTIDES = ["A", "C", "G", "T"]
RNA_NUCLEOTIDES = ["A", "C", "G", "U"]


# compter les nucléotides dans l'ADN
def CountNucletoides(dna_string):
    return [ (c, dna_string.count(c)) for c in DNA_NUCLEOTIDES]



def DNAtoRNA(dna_string):   
    '''
        simplement replacer les T (thymine) pas des U (uracile)
    '''
    return dns_string.replace("T", "U")


def ReverseComplement(dna_string):
    '''
        inverser, puis swapper T<->A et C<->G
    '''
    return dna_string[::-1].replace("T", "Z").replace("A", "T").replace("Z", "A").replace("C", "Z").replace("G", "C").replace("Z", "G")
