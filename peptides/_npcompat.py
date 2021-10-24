import functools
import operator


def prod(a):
    """Return the product of an iterable of numbers
    """
    return functools.reduce(operator.prod, a)
