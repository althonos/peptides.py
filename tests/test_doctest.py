# coding: utf-8
"""Test doctest contained tests in every file of the module.
"""

import configparser
import doctest
import importlib
import os
import pkgutil
import re
import sys
import shutil
import types
import warnings
from unittest import mock

try:
    import numpy
except ImportError:
    numpy = None

import peptides


def _load_tests_from_module(tests, module, globs, setUp=None, tearDown=None):
    """Load tests from module, iterating through submodules.
    """
    for attr in (getattr(module, x) for x in dir(module) if not x.startswith("_")):
        if isinstance(attr, types.ModuleType):
            suite = doctest.DocTestSuite(
                attr,
                globs,
                setUp=setUp,
                tearDown=tearDown,
                optionflags=+doctest.ELLIPSIS,
            )
            tests.addTests(suite)
    return tests


def load_tests(loader, tests, ignore):
    """`load_test` function used by unittest to find the doctests.
    """

    def setUp(self):
        warnings.simplefilter("ignore")
        if numpy is not None:
            numpy.set_printoptions(legacy='1.25')

    def tearDown(self):
        warnings.simplefilter(warnings.defaultaction)

    def setUp2(self):
        setUp(self)
        self._patch = mock.patch("peptides.numpy", new=None)
        self._patch.__enter__()

    def tearDown2(self):
        tearDown(self)
        self._patch.__exit__(None, None, None)

    # doctests are not compatible with `green`, so we may want to bail out
    # early if `green` is running the tests
    if sys.argv[0].endswith("green"):
        return tests

    # recursively traverse all library submodules and load tests from them
    packages = [None, peptides]

    for pkg in iter(packages.pop, None):
        # import the base module and add it to the tests
        globs = dict(peptides=peptides, **pkg.__dict__)
        # run the tests using either `
        tests.addTests(
            doctest.DocTestSuite(
                pkg,
                globs=globs,
                setUp=setUp,
                tearDown=tearDown,
                optionflags=+doctest.ELLIPSIS,
            )
        )
        # if we ran the doctests with NumPy, run the doctests a second
        # time, but this time using the fallback implementation
        if numpy is not None:
            tests.addTests(
                doctest.DocTestSuite(
                    pkg,
                    globs=globs,
                    setUp=setUp2,
                    tearDown=tearDown2,
                    optionflags=+doctest.ELLIPSIS,
                )
            )
        # find submodules
        for (_, subpkgname, subispkg) in pkgutil.walk_packages(pkg.__path__):
            # do not import __main__ module to avoid side effects!
            if subpkgname == "__main__":
                continue
            # if the submodule is a package, we need to process its submodules
            # as well, so we add it to the package queue
            if subispkg:
                module = importlib.import_module(".".join([pkg.__name__, subpkgname]))
                packages.append(module)

    return tests
