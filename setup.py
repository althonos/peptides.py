#!/usr/bin/env python
# coding: utf-8

import ast
import csv
import configparser
import glob
import math
import os

import setuptools
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.build_py import build_py as _build_py

try:
    import astor
except ImportError as err:
    astor = err


class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class codegen(setuptools.Command):
    """A custom command to build Python sources from CSV data files.
    """

    description = "build source code from CSV data tables"
    user_options = [
        ("force", "f", "force rebuilding the files even if they are not outdated"),
        ("inplace", "i", "copy the build files to the source directory when done"),
        ("tables",  "t", "the path to the data folder containing the CSV tables"),
        ("datasets",  "d", "the path to the data folder containing the datasets"),
    ]

    def initialize_options(self):
        self.force = False
        self.inplace = False
        self.tables = None
        self.datasets = None

    def finalize_options(self):
        _build_py = self.get_finalized_command("build_py")
        self.build_lib = _build_py.build_lib
        if self.tables is None:
            self.tables = os.path.relpath(os.path.join(__file__, "..", "peptides", "tables"))
        if self.datasets is None:
            self.datasets = os.path.relpath(os.path.join(__file__, "..", "peptides", "datasets"))

    # ----------------------------

    def _load_csv(self, filename):
        self.announce(f"loading {filename!r}", level=1)
        with open(filename, "r") as f:
            return {row[0]:float(row[1]) for row in csv.reader(f)}

    def _load_matrix(self, filename):
        self.announce(f"loading {filename!r}", level=1)
        with open(filename, "r") as f:
            return [
                [float(x) for x in line.split()]
                for line in f
                if not line.startswith("#")
            ]

    def _load_tables(self, path):
        self.announce(f"loading data tables from {path!r}", level=2)
        tables = {}
        for group in os.listdir(path):
            if group == "__pycache__":
                continue
            if not os.path.isdir(os.path.join(path, group)):
                continue
            tables[group] = {}
            for filename in glob.iglob(os.path.join(path, group, "*.csv")):
                member, _ = os.path.splitext(os.path.basename(filename))
                tables[group][member] = self._load_csv(filename)
        return tables

    def _load_dataset(self, path):
        dataset = self._load_tables(path)
        for group in os.listdir(path):
            for filename in glob.iglob(os.path.join(path, group, "eigen.txt")):
                eigenmatrix = self._load_matrix(filename)
                dataset[group]["eigenvalues"] = [m[0] for m in eigenmatrix]
                dataset[group]["eigenvectors"] = [m[1:] for m in eigenmatrix]
        return dataset

    def _load_datasets(self, path):
        self.announce(f"loading datasets from {path!r}", level=2)
        datasets = {}
        for dataset in os.listdir(path):
            if dataset == "__pycache__":
                continue
            if not os.path.isdir(os.path.join(path, dataset)):
                continue
            datasets[dataset] = self._load_dataset(os.path.join(self.datasets, dataset))
        return datasets


    # ----------------------------

    def _make_literal(self, value):
        if isinstance(value, (int, float, str)):
            return ast.Constant(value)
        elif isinstance(value, dict):
            return ast.Dict(
                keys=[self._make_literal(x) for x in value.keys()],
                values=[self._make_literal(x) for x in value.values()]
            )
        elif isinstance(value, list):
            return ast.List(
                elts=[self._make_literal(x) for x in value]
            )
        else:
            raise TypeError(f"cannot make literal for {value!r}")

    def _generate_tables_module(self, tables):
        n = sum(map(len, tables.values()))
        self.announce(f"building Python AST for {n!r} data tables", level=2)
        body = []
        for name, table in tables.items():
            stub = name.upper().replace(".", "_")
            body.append(ast.Assign(
                targets=[ast.Name(id=stub, ctx=ast.Store())],
                value=self._make_literal(table)
            ))
        return ast.Module(body=body)

    def _generate_datasets_module(self, datasets):
        body = []
        for name, dataset in datasets.items():
            n = sum(map(len, dataset.values()))
            self.announce(f"building Python AST for the {n!r} data tables of the {name!r} dataset", level=2)
            stub = name.upper().replace(".", "_")
            body.append(ast.Assign(
                targets=[ast.Name(id=stub, ctx=ast.Store())],
                value=self._make_literal(dataset)
            ))
        return ast.Module(body=body)

    def _write_module(self, module, filename):
        self.announce(f"writing Python source to {filename!r}", level=2)
        with open(filename, "w") as f:
            f.write("# this file was automatically generated\n")
            f.write("# by the `python setup.py codegen` command\n")
            f.write("# DO NOT EDIT MANUALLY!\n")
            f.write(astor.to_source(module))
            f.write("\n")
        if self.inplace:
            localpath = os.path.relpath(filename, start=self.build_lib)
            self.copy_file(filename, localpath)

    def run(self):
        # check astor is available
        if isinstance(astor, ImportError):
            raise RuntimeError("`astor` package is required for the `codegen` command") from astor

        # generate a Python file containing the constants tables
        tables = self._load_tables(self.tables)
        tables_module = self._generate_tables_module(tables)
        table_file = os.path.join(self.build_lib, "peptides", "tables", "__init__.py")
        self.mkpath(os.path.dirname(table_file))
        self._write_module(tables_module, table_file)

        # generate a Python file containing some datasets with nested / diverse data
        datasets = self._load_datasets(self.datasets)
        datasets_module = self._generate_datasets_module(datasets)
        datasets_file = os.path.join(self.build_lib, "peptides", "datasets", "__init__.py")
        self.mkpath(os.path.dirname(datasets_file))
        self._write_module(datasets_module, datasets_file)


class build_py(_build_py):
    """A hacked `build_py` command that will also run `codegen`.
    """

    def run(self):
        # generate tables if needed
        if not self.distribution.have_run.get("codegen", False):
            _codegen = self.get_finalized_command("codegen")
            _codegen.force = self.force
            _codegen.run()
        # build rest as normal
        _build_py.run(self)


if __name__ == "__main__":
    setuptools.setup(cmdclass=dict(build_py=build_py, codegen=codegen))
