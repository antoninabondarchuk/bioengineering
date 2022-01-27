"""This module is dedicated to the very basics of
finding origin of replication in DNA.

Content:
    1. find_minimum_skew() function.

The supplementary materials can be found on
the correspondent Coursera/Stepik course here:
https://www.coursera.org/learn/dna-analysis

Week 2.
https://www.coursera.org/learn/dna-analysis/home/week/2

Author:
Antonina Bondarchuk
antonina.bondarchuk1@gmail.com
2022

"""


def find_minimum_skew(genome):
    """Iterates through genome and searches for the
    minimums while counting the 'level', where decrements
    if faces C and increments if faces G.

    Finding a minimums can give us the position of where
    origin of replication is. Quantity of G minus C
    tarts to increase after reaching ori.

    Notes:
        Iteration assumes that 0 position has level 0,
        therefore it starts from index 1 which is
        actually the 1st char in genome.

    Examples:
        >>> find_minimum_skew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
        [11, 24]

    Args:
        genome (str): A DNA sequence to process.

    Returns:
        list[int]. Indexes (counting from 1) of
        minimums in genome, possible ori.

    """
    curr_level = 0
    start_idx = None
    min_levels_idxs_pairs = [(curr_level, start_idx), ]
    for i, char in enumerate(genome, start=1):
        if char == 'C':
            curr_level -= 1

            curr_min_level = min_levels_idxs_pairs[0][0]
            if curr_level == curr_min_level:
                min_levels_idxs_pairs.append((curr_level, i))
            elif curr_level < curr_min_level:
                min_levels_idxs_pairs.clear()
                min_levels_idxs_pairs.append((curr_level, i))

        elif char == 'G':
            curr_level += 1

    min_level_indexes = [idx for lvl, idx in min_levels_idxs_pairs]
    return min_level_indexes
