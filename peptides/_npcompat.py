"""Compatibility layer to use `array.array` when `numpy` is not available.
"""

import array
import functools
import operator


class array(array.array):

    def __new__(cls, values, dtype='f'):
        return super().__new__(cls, dtype, values)


def take(a, indices):
    """Take elements from an array.
    """
    return array([a[i] for i in indices])


def prod(a):
    """Return the product of an iterable of numbers.
    """
    return functools.reduce(operator.prod, a)


def zeros(shape, dtype='f'):
    return array([0 for _ in range(shape)], dtype=dtype)
