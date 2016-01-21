import random

from itertools import repeat, chain

import numpy as np
from numpy.random import beta
from numpy.random import exponential


def phred(p):
    return -10 * np.log10(p)


def sample_error_ps(length, errors):
    """Probability of error for each base."""
    mean = errors / length
    a = 1
    b = (1 - mean) / mean
    return beta(a, b, size=length)
    

def sample_expected_errors(mean, n):
    """Expected number of errors for a read"""
    return exponential(mean, size=n)


def sample_template(length):
    return ''.join(random.choice('CGTA') for _ in range(length))


def mutate_base(base):
    return random.choice(list(set('ACGT') - set(base)))


def mutate_template(template, identity):
    max_errors = np.floor(len(template) * (1 - identity) / 2)
    if max_errors < 1:
        raise Exception('identity is too small: {}'.format(identity))
    n_errs = random.randint(1, max_errors)
    positions = list(sorted(random.sample(range(len(template)), n_errs)))
    mutations = list(mutate_base(template[p]) for p in positions)
    result = list(template)
    for i, c in zip(positions, mutations):
        result[i] = c
    return ''.join(result)


def coinflip(p):
    return random.random() < p


def read_base(base, p):
    if coinflip(p):
        return mutate_base(base)
    return base


def generate_read(template, ps):
    return ''.join(list(read_base(c, p) for c, p in zip(template, ps)))


def repeat_and_chain(items, ns):
    """
    
    >>> list(repeat_and_chain([3, 5], [2, 2]))
    [3, 3, 5, 5]

    """
    return chain.from_iterable(repeat(item, n) for item, n in zip(items, ns))


def sample(identity, ns, length=1000, error_mean=5, shuffle=True):
    """Generate similar templates, and sample from them.

    Paramters
    ---------
    identity : float
        In range [0, 1]. Minimum pairwise identity for all templates.

    ns : list of integers
        Number of sequences to sample from each template.

    length : integer
        Length of sequences

    error_mean : float
        Mean number of expected errors per sequence.

    Examples
    --------
    sample(0.99, [50, 50, 50], length=1000, error_mean=10)

    """
    base_template = sample_template(length)
    templates = list(mutate_template(base_template, identity) for _ in range(len(ns)))
    errors = sample_expected_errors(error_mean, sum(ns))
    ps = np.vstack(list(sample_error_ps(length, e) for e in errors))
    reads = list(generate_read(template, qv) for template, qv in zip(repeat_and_chain(templates, ns), ps))
    template_indices = list(repeat_and_chain(range(len(templates)), ns))
    if shuffle:
        indices = list(range(sum(ns)))
        random.shuffle(indices)
        ps = ps[indices]
        reads = list(reads[i] for i in indices)
        template_indices = list(template_indices[i] for i in indices)
    return templates, phred(ps), reads, template_indices
