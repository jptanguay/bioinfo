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


def FB(n,k,i=0,v=1):
    v += k
    if i >= n:
        return v    
    return fu(n,k,i+1,v)
    
    
# pas encore au point 
import Levenshtein 
def mt(ge,ds,k):
    
    """
        pip3 install Levenshtein
        
        
        Sample Dataset
        2
        ACGTAG
        ACGGATCGGCATCGT
        Sample Output
        1 4
        1 5
        1 6
    
    """
    for i in range( len(ge)-len(ds) ):
        for j in range(k+1):
            #print(ds[0:len(ds)-j])
            d = Levenshtein.distance(ds[0:len(ds)-j], ge[i:len(ds)])
            #print("d", d)
            if d <= k:
                print(i+1, len(ds)-j)
                
                
def Hamming(s,t):
    #s = "GAGCCTACTAACGGGAT"
    #t = "CATCGTAATGACGGCCT"

    d = sum([a!=b for (a,b) in zip(s,t)])
    return d


import itertools
def permu(n):
    r = [i+1 for i in range(n)]
    perms = list(itertools.permutations(r))
    #print(perms)
    print(len(perms))
    for p in perms:
        s = " ".join(str(item) for item in p)
        print(s)
        