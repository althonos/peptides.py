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

    def tearDown(self):
        warnings.simplefilter(warnings.defaultaction)

    # doctests are not compatible with `green`, so we may want to bail out
    # early if `green` is running the tests
    if sys.argv[0].endswith("green"):
        return tests

    # doctests require `numpy` to run, which may not be available because
    # it is a pain to get to work out-of-the-box on OSX
    # if numpy is None:
    #     return tests

    # recursively traverse all library submodules and load tests from them
    packages = [None, peptides]

    for pkg in iter(packages.pop, None):
        # import the base module and add it to the tests
        globs = dict(numpy=numpy, peptides=peptides, **pkg.__dict__)
        tests.addTests(
            doctest.DocTestSuite(
                pkg,
                globs=globs,
                setUp=setUp,
                tearDown=tearDown,
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
