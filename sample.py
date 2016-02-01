import random

from itertools import repeat, chain

import numpy as np
from numpy.random import poisson


def phred(p):
    return -10 * np.log10(p)


def sample_template(length):
    return ''.join(random.choice('CGTA') for _ in range(length))


def mutate_base(base):
    return random.choice(list(set('ACGT') - set(base)))


def coinflip(p):
    return random.random() < p


def sample_from_template(template, point_rate, insertion_rate, deletion_rate):
    result = []
    for base in template:
        while coinflip(insertion_rate):
            result.append(random.choice('ACGT'))
        if coinflip(deletion_rate):
            continue
        if coinflip(point_rate):
            result.append(mutate_base(base))
        else:
            result.append(base)
    while coinflip(insertion_rate):
        result.append(random.choice('ACGT'))
    return ''.join(result)
        

def repeat_and_chain(items, ns):
    """
    
    >>> list(repeat_and_chain([3, 5], [2, 2]))
    [3, 3, 5, 5]

    """
    return chain.from_iterable(repeat(item, n) for item, n in zip(items, ns))


def sample(n, length, mean_error_rate, insertion_rate, deletion_rate):
    """Generate a template, and sample from it.

    Paramters
    ---------
    s : list of integers
        Number of sequences to sample from each template.

    length : integer
        Length of sequences

    mean_error_rate : float
        Mean number of expected point errors per sequence.

    """
    reads = []
    ps = []
    template = sample_template(length)
    for i in range(n):
        point_rate = poisson(mean_error_rate)
        reads.append(sample_from_template(template,
                                          point_rate / length,
                                          insertion_rate / length,
                                          deletion_rate / length))
        ps.append(list(repeat(point_rate / length, length)))
    return template, reads, ps
