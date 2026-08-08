#!/usr/bin/env python

# Author: <Imre Rist> <imre@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.6"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = "ATCGNatcgn"
RNA_bases = "AUCGNaucgn"
DNA_complement = str.maketrans('ATCGatcg', 'TAGCTAGC')
RNA_complement = str.maketrans('AUCGaucg', 'UAGCUAGC')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Returns the average quality score of a given phred score sequence string"""
    score_sums = 0
    for char in phred_score:
        score_sums = score_sums + convert_phred(char)
    return score_sums / len(phred_score)

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return seq.strip(RNA_bases if RNAflag else DNA_bases) == ""

def gc_content(seq: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq, True) or validate_base_seq(seq, False), 'Not a valid DNA or RNA sequence'
    return (seq.upper().count("G") + seq.upper().count("C")) / len(seq)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    if len(lst) % 2 == 0: # Even number of elements, need to average 2 middle values
        mid_ciel = len(lst) // 2
        return (lst[mid_ciel] + lst[mid_ciel - 1]) / 2
    else:
        return lst[len(lst) // 2]

def oneline_fasta(in_file: str, out_file: str) -> None:
    '''Takes an input fasta file, writes the same information to the output file with each sequence only occupying one line'''
    with open(out_file, 'w') as out_fh:
        with open(in_file, 'r') as in_fh:
            first_line = True
            for line in in_fh:
                if first_line:
                    out_fh.write(line)
                    first_line = False
                elif line[0] == '>':
                    out_fh.write(f'\n{line}')
                else:
                    out_fh.write(line.strip())

def reverse_complement(seq: str, RNAflag: bool = False) -> str:
    '''Takes a sequence string, returns the reverse complement. Uses U instead of T if RNAflag is True.
    Note: To save time, this function does not validate the sequence.
    Ensure the input sequence is a valid DNA or RNA sequence beforehand.'''
    return seq.translate(RNA_complement if RNAflag else DNA_complement)[::-1]


if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)

    # Test convert_phred():
    assert convert_phred("I") == 40, "Wrong phred score for 'I'"
    assert convert_phred("C") == 34, "Wrong phred score for 'C'"
    assert convert_phred("2") == 17, "Wrong phred score for '2'"
    assert convert_phred("@") == 31, "Wrong phred score for '@'"
    assert convert_phred("$") == 3, "Wrong phred score for '$'"
    print("convert_phred() function is working")

    # Test qual_score():
    assert qual_score("IIII#I#####CCCCC") == 23.875
    assert qual_score('#####################') == 2
    print("qual_score() function is working")

    # Test validate_base_seq():
    assert validate_base_seq('ATgatcagAGTcaA', False)
    assert not validate_base_seq('ATgatcagAGTcaA', True)
    assert validate_base_seq('AuAUCgaGACaAUGg', True)
    assert not validate_base_seq('AuAUCgaGACaAUGg', False)
    assert not validate_base_seq('ACGATCLa', False)
    assert not validate_base_seq('ACGATCLa', True)
    print("validate_base_seq() function is working")

    # Test gc_content():
    assert gc_content("ATTAGCTAGCTAGGA") == 0.4, "Wrong GC content for ATTAGCTAGCTAGGA"
    assert gc_content("tagctagagc") == 0.5, "Wrong GC content for tagctagagc"
    assert gc_content("TAgcatGCTAtgacGTCAGT") == 0.45, "Wrong GC content for TAgcatGCTAtgacGTCAGT"
    assert gc_content("GGGGAGCcagGT") == 0.75, "Wrong GC content for GGGGAGCcagGT"
    assert gc_content("ATTATATAtATATATAaaaTA") == 0.0, "Wrong GC content for ATTATATAtATATATAaaaTA"
    assert gc_content("GCGCGCccGCgCGggCGCgg") == 1.0, "Wrong GC content for GCGCGCccGCgCGggCGCgg"
    print("gc_content() function is working")
    
    # Test calc_median():
    assert calc_median([1,2,3,4,5]) == 3
    assert calc_median([2,6,14,33,51,105]) == 23.5
    print("calc_median() function is working")

    # Test reverse_complement():
    assert reverse_complement('ATaGCAtaCGG') == 'CCGTATGCTAT'
    assert reverse_complement('NTAcgNT') == 'ANCGTAN'
    assert reverse_complement('NNTAGANC') == 'GNTCTANN'
    print("reverse_complement() function is working")
