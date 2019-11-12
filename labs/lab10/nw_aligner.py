#!/usr/bin/python3

"""

Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

"""
import sys
import os

if sys.version_info[0] != 3:
    print('Please use python3 to run this script.')
    sys.exit(1)

from Bio import SeqIO
import numpy as np

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        """
        Input: (String) A path to a scoring matrix file.
        Output: (Dictionary) A nested dictionary where keys are strings
                and elements are scores as integers.
    
        Example:
    
        >>> matrix, gap_penalty = NWAligner.load_score_matrix('/home/bioe131/BLOSUM62')
        >>> matrix['A']['A']
        4
        >>> matrix['W']['W']
        11
        >>> gap_penalty
        -4

        """

        score_matrix = {}
        gap_penalty = None
        keys = []
        non_comment_line_count = 0
        with open(fname) as fp:
            for line_num, line in enumerate(fp):
                # ignore comments in matrix file
                if line.startswith("#"):
                    continue
                
                ### TODO ### - DONE
                # Parse matrix file line-by-line and load into nested dictionaries.
                # 
                # Last line of matrix contains the gap penalty which must be pulled
                # out and returned.
                if non_comment_line_count == 0: # we should create our keys
                    keys = line.strip().split(" ")
                elif non_comment_line_count == len(keys) + 1:
                    gap_penalty = int(line.strip())
                    score_matrix["-"] = {}
                    for i in range(len(keys)):
                        score_matrix["-"][keys[i]] = gap_penalty
                        score_matrix[keys[i]]["-"] = gap_penalty
                else:
                    score_matrix[keys[non_comment_line_count - 1]] = {}
                    line = line.strip().split(" ")
                    for i in range(len(line)):
                        score_matrix[keys[non_comment_line_count - 1]][keys[i]] = int(line[i])
                
                non_comment_line_count += 1
                
        return score_matrix, gap_penalty

    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """

        seqs = []

        ### TODO ### - DONE
        # Load FASTA file and return list of sequences.
        # Throw an error if there are more than two sequences in the file.
        
        records = list(SeqIO.parse(fname, "fasta"))
        if len(records) != 2:
            raise Exception('FASTA file does not contain exactly two sequences.')
        seqs = [str(record.seq) for record in records]
        
        return seqs

    def align(self, seq_x, seq_y, print_matrix = False):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]

        ### TODO ### - DONE
        # Fill the top row of the matrix with scores for gaps.
        # Fill the first column of the matrix with scores for gaps.
        
        gap_penalty_increasing = 0
        for i in range(0, max(len(seq_x), len(seq_y)) + 1):
            if i <= len(seq_x):
                matrix[i][0] = gap_penalty_increasing
            if i <= len(seq_y):
                matrix[0][i] = gap_penalty_increasing
            gap_penalty_increasing += self.gap_penalty

        ###
        ### RECURSION
        ###

        # fill the dynamic programming and pointer matrices
        for x in range(1, len(seq_x) + 1):
            for y in range(1, len(seq_y) + 1):
                match_score = self.score_matrix[seq_x[x-1]][seq_y[y-1]]
                
                ### TODO ### - DONE
                # Take the maximum score of three possibilities:
                #   1) The element in the matrix diagonal from this one
                #      plus the score of an exact match
                #   2) The element to the left plus a gap penalty
                #   3) The element above plus a gap penalty
                # ... and set the current element (matrix[x][y]) equal to that
                #
                # Keep track of which of these choices you made by setting
                #   the same element (i.e., pointers[x][y]) to some value that
                #   has meaning to you.
                
                values = [
                    matrix[x-1][y-1] + match_score,
                    matrix[x-1][y]   + self.gap_penalty,
                    matrix[x][y-1]   + self.gap_penalty
                ]
                pointers[x][y] = np.argmax(values) + 1
                matrix[x][y] = np.max(values)

        # print the dynamic programming matrix
        if print_matrix:
            for x in range(len(seq_x) + 1):
                print("\t".join(map(lambda i: str(int(i)), matrix[x])))
                
            for x in range(len(seq_x) + 1):
                print("\t".join(map(lambda i: str(int(i)), pointers[x]))) 

        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x), len(seq_y)

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []

        while x > 0 or y > 0:
            move = pointers[x][y]

            ### TODO ###
            # Follow pointers back through the matrix to the origin.
            # Depending on which "move" you made at each element in the
            #   matrix, you'll either align seq_x to seq_y, seq_x to a gap, or
            #   seq_y to a gap.
            
            if move == 0:
                # at the beginning
                if x == 0:
                    align_x.append("-")
                    align_y.append(seq_y[y-1])
                    y -= 1
                else: # y == 0
                    align_x.append(seq_x[x-1])
                    align_y.append("-")
                    x -= 1
            elif move == 1:
                # from the diagonal
                align_x.append(seq_x[x-1])
                align_y.append(seq_y[y-1])
                x -= 1
                y -= 1
            elif move == 2:
                # the character in sequence N was aligned to a gap, equivalent to move up
                align_x.append(seq_x[x-1])
                align_y.append("-")
                x -= 1
            elif move == 3:
                # the character in sequence M was aligned to a gap, equivalent to move left
                align_x.append("-")
                align_y.append(seq_y[y-1])
                y -= 1
        self.matrix = matrix
        self.pointers = pointers
        # flip the alignments, as they're reversed
        self.seq1 = "".join(align_x[::-1])
        self.seq2 = "".join(align_y[::-1])
        return (self.seq1, self.seq2)
    
    def score_seqs(self, seq1, seq2):
        """
        Input: (Strings) Two sequences that have been aligned (seq_x and seq_y).
        Output: (Tuple of strings) The score, given the aligner's score_matrix.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---', 'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')
        >>> aligner.score_seqs(seqs[0], seqs[1])
        16
        """
        score = 0
        for i in range(len(seq1)):
            score += self.score_matrix[seq1[i]][seq2[i]]
        return score
        

###                                      ###
### NO NEED TO EDIT CODE BELOW THIS LINE ###
###                                      ###

if __name__ == '__main__':
    def usage():
        print('usage: %s matrixfilename stringfilename')
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print('Can not open %s' % (fname,))
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))
