"""This module is dedicated to the very basics of
finding origin of replication in DNA.

Content:
    1. find_minimum_skew() function.
    2. hamming_distance() function.

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


def hamming_distance(str_a, str_b):
    """Calculates the number of characters which
    are pairwise different in the given string of
    equal length.

    Examples:
        >>> hamming_distance('GGGCCGTTGGT', 'GGACCGTTGAC')
        3

    See Also:
        https://en.wikipedia.org/wiki/Hamming_distance

    Args:
        str_a (str): A string to compare with str_b.
        str_b (str): A string to compare with str_a.

    Returns:
        distance (int). Quantity of characters differ.

    """
    distance = 0
    for i in range(len(str_a)):
        if str_a[i] != str_b[i]:
            distance += 1
    return distance


def approximate_pattern_matching_starts(pattern, genome, max_diff):
    """Finds all approximate occurrences of a pattern in a string
    while iterating through genome - pattern length.

    References:
        hamming_distance()

    Examples:
        >>> approximate_pattern_matching_starts('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTC', 3)
        [6, 7, 26, 27]

    Args:
        pattern (str): Substring to match.
        genome (str): A DNA sequence to process.
        max_diff (int): Maximum number of mismatched
        characters between pattern and substring of genome.

    Returns:
        matching_start_positions (list[int]). List of start
        indexes of substrings that fully or slightly
        (according to max_diff) matches in ascending order.

    """
    pattern_len = len(pattern)
    genome_len_to_iterate = len(genome) - pattern_len + 1
    matching_start_positions = list()
    for i in range(genome_len_to_iterate):
        if hamming_distance(genome[i:i+pattern_len], pattern) <= max_diff:
            matching_start_positions.append(i)
    return matching_start_positions
