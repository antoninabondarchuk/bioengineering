"""This module is dedicated to the very basics of
finding hidden messages in DNA.

Content:
    1. get_count_pattern() function.
    2. get_frequency_table() function.

The supplementary materials can be found on
the correspondent Coursera/Stepik course here:
https://www.coursera.org/learn/dna-analysis

Week 1.
https://www.coursera.org/learn/dna-analysis/home/week/1

Author:
Antonina Bondarchuk
antonina.bondarchuk1@gmail.com

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
    as slice ([:]) function excludes the upper boundary.

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
