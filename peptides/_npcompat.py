"""Compatibility layer to use `array.array` when `numpy` is not available.
"""

import array as _array
import builtins
import functools
import operator

def array(values, dtype="f"):
    _values = _array.array(dtype, values)
    return lazyarray(_values, dtype=dtype, length=len(_values))



class lazyarray(object):

    def __init__(self, values, dtype, length):
        self.__values = values
        self.dtype = dtype
        self.length = length

    def __iter__(self):
        return iter(self.__values)

    def __len__(self):
        return self.length

    def __getitem__(self, i):
        if isinstance(i, slice):
            s = self.__values.__getitem__(i)
            return lazyarray(s, self.dtype, len(s))
        return self.__values.__getitem__(i)

    def __setitem__(self, i, val):
        self.__values.__setitem__(i, val)

    def __add__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((x + other for x in self), self.dtype, self.length)
            # return lazyarray((x + other for x in self), self.dtype, self.length)
        elif not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x1 + x2 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __radd__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((other + x for x in self), self.dtype, self.length)
        elif not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x2 + x1 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((x - other for x in self), self.dtype, self.length)
        elif not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x1 - x2 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __rsub__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((other - x for x in self), self.dtype, self.length)
        elif not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x2 - x1 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((x*other for x in self), self.dtype, self.length)
        elif not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x1 * x2 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((other*x for x in self), self.dtype, self.length)
        elif not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x2 * x1 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((x/other for x in self), self.dtype, self.length)
        elif not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x1 / x2 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((other/x for x in self), self.dtype, self.length)
        elif not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x2 / x1 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __pow__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((x**other for x in self), self.dtype, self.length)
        if not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x1 ** x2 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def __rpow__(self, other):
        if isinstance(other, (int, float)):
            return lazyarray((other**x for x in self), self.dtype, self.length)
        if not isinstance(other, lazyarray):
            return NotImplemented
        if len(self) != len(other):
            raise ValueError("Cannot pairwise multiply arrays of different lengths")
        return lazyarray(
            (x2 ** x1 for x1, x2 in zip(self, other)),
            self.dtype,
            self.length
        )

    def sum(self):
        return builtins.sum(self.__values)

    def take(self, indices):
        return lazyarray(
            map(self.__values.__getitem__, indices),
            self.dtype,
            len(indices)
        )

    def sum(self):
        return builtins.sum(self.__values)


def take(a, indices):
    """Take elements from an array.
    """
    return a.take(indices)


def prod(a):
    """Return the product of an iterable of numbers.
    """
    return functools.reduce(operator.prod, a)


def zeros(shape, dtype='f'):
    return lazyarray([0 for _ in range(shape)], dtype=dtype, length=shape)
