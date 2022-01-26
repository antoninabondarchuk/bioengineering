"""This module is dedicated to the very basics of
finding hidden messages in DNA.

Content:
    1. get_count_pattern() function.

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
