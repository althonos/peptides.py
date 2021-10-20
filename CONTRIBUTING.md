# Contributing to `peptides.py`

For bug fixes or new features, please file an issue before submitting a
pull request. If the change isn't trivial, it may be best to wait for
feedback.

## Running tests

Tests are written as usual Python unit tests with the `unittest` module of
the standard library. Running them requires the data files to be built
locally:

```console
$ python setup.py codegen --inplace
$ python -m unittest discover -vv
```

## Coding guidelines

This project targets Python 3.6 or later.

