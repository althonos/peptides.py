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

### Docstrings

The docstring lines should not be longer than 76 characters (which allows rendering without soft-wrap of the entire module in a 80x24 terminal window).  
Docstrings should be written in Google format.
