"""Compatibility layer to use `array.array` when `numpy` is not available.
"""

import array
import functools
import operator


class array(array.array):

    def __new__(cls, values, dtype='f'):
        return super().__new__(cls, dtype, values)

    def __add__(self, other):
        if isinstance(other, (int, float)):
            return array((x + other for x in self), dtype=self.typecode)
        elif not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x1 + x2 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __radd__(self, other):
        if isinstance(other, (int, float)):
            return array((other + x for x in self), dtype=self.typecode)
        elif not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x2 + x1 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return array((x - other for x in self), dtype=self.typecode)
        elif not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x1 - x2 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __rsub__(self, other):
        if isinstance(other, (int, float)):
            return array((other - x for x in self), dtype=self.typecode)
        elif not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x2 - x1 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return array((x*other for x in self), dtype=self.typecode)
        elif not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x1 * x2 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return array((other*x for x in self), dtype=self.typecode)
        elif not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x2 * x1 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return array((x/other for x in self), dtype=self.typecode)
        elif not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x1 / x2 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            return array((other/x for x in self), dtype=self.typecode)
        elif not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x2 / x1 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __pow__(self, other):
        if isinstance(other, (int, float)):
            return array((x**other for x in self), dtype=self.typecode)
        if not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x1 ** x2 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )

    def __rpow__(self, other):
        if isinstance(other, (int, float)):
            return array((other**x for x in self), dtype=self.typecode)
        if not isinstance(other, array):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return array(
            (x2 ** x1 for x1, x2 in zip(self, other)),
            dtype=self.typecode
        )


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
