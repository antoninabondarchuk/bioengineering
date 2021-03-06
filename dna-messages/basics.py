"""This module is dedicated to the very basics of
finding hidden messages in DNA.

Content:
    1. get_count_pattern() function.
    2. get_frequency_table() function.
    3. get_most_frequent_patterns() function.
    4. get_complementary_sequence() function.
    5. get_pattern_start_positions() function.

The supplementary materials can be found on
the correspondent Coursera/Stepik course here:
https://www.coursera.org/learn/dna-analysis

Week 1.
https://www.coursera.org/learn/dna-analysis/home/week/1

Author:
Antonina Bondarchuk
antonina.bondarchuk1@gmail.com
2022

"""


def get_count_pattern(text, pattern):
    """Counts how many occurrences of pattern in text.
    Counts even overlapping occurrences.

    Iterates from index 0 to length of text
    minus length of the pattern and plus one
    to include the last possible match
    as range function excludes the upper boundary.

    Examples:
        >>> get_count_pattern('GACCATCAAAACTGATAAACTACTTAAAAATCAGTAAA', 'AAA')
        7

    Args:
        text (str): DNA sequence.
        pattern (str): part of the DNA sequence.

    Returns:
        count (int).
        The quantity of times pattern matched text part.

    """
    pattern_len = len(pattern)
    text_len_to_iterate = len(text) - pattern_len
    count = 0
    for i in range(text_len_to_iterate + 1):
        if text[i:i + pattern_len] == pattern:
            count += 1
    return count


def get_frequency_table(text, k):
    """Calculates the quantities of occurrences
    of all possible patterns of length k inside the text.

    Iterates from index 0 to length of text
    minus length of the pattern and plus one
    to include the last possible match
    as range function excludes the upper boundary.

    Examples:
        >>> get_frequency_table('CAAAAACTCAAA', 3)
        {'CAA': 2, 'AAA': 4, 'AAC': 1, 'ACT': 1, 'CTC': 1, 'TCA': 1}

    Args:
        text (str): The DNA sequence to analyse.
        k (int): Number to find k-mers (patterns of length k).

    Returns:
        frequencies (dict).
        Mapping with patterns as keys and their counts as values.

    """
    frequencies = dict()
    text_len_to_iterate = len(text) - k + 1
    for i in range(text_len_to_iterate):
        pattern = text[i:i + k]
        if pattern not in frequencies:
            frequencies[pattern] = 1
        else:
            frequencies[pattern] += 1
    return frequencies


def get_most_frequent_patterns(text, k):
    """Calculates frequency table for each pattern
    of length k in text and takes those,
    which numbers of occurrencies in the text
    are the maximum ones.

    Examples:
        >>> get_most_frequent_patterns('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4)
        ['GCAT', 'CATG']

    Args:
        text (str): The DNA sequence.
        k (int): Number to find k-mers (patterns of length k).

    Returns:
        most_freq_patterns (list[str]).
        List of the most frequently occurred patterns
        in text of length k.

    """
    frequencies = get_frequency_table(text, k)
    max_occurred_num = max(frequencies.values())
    most_freq_patterns = list()
    for pattern, occurrences in frequencies.items():
        if occurrences == max_occurred_num:
            most_freq_patterns.append(pattern)
    return most_freq_patterns


def get_complementary_sequence(pattern):
    """Reverses pattern and places a complementary
    nucleotide from the complementary_nucleotides dictionary
    to the result.

    Examples:
        >>> get_complementary_sequence('AAAACCCGGT')
        'ACCGGGTTTT'

    Args:
        pattern (str): The DNA sequence to process.

    Returns:
        complement (str).
        Complementary DNA sequence to the given pattern.

    """
    complementary_nucleotides = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = ''
    for char in pattern[::-1]:
        complement += complementary_nucleotides[char]
    return complement


def get_pattern_start_positions(pattern, text):
    """Iterates from index 0 to length of text
    minus length of the pattern and plus one
    to include the last possible match
    as range function excludes the upper boundary.

    Finds all matching substrings in text, and writes
    their start positions.

    Examples:
        >>> get_pattern_start_positions('ATAT', 'GATATATGCATATACTTATAT')
        [1, 3, 9, 17]

    Args:
        pattern (str): Substring to find in text.
        text (str): The DNA sequence to process.

    Returns:
        pattern_start_positions (list[int]).
        List of starting positions of pattern
        in the text in ascending order.

    """
    len_pattern = len(pattern)
    pattern_start_positions = list()
    text_len_to_iterate = len(text) - len_pattern + 1
    for i in range(text_len_to_iterate):
        if text[i:i + len_pattern] == pattern:
            pattern_start_positions.append(i)
    return pattern_start_positions
