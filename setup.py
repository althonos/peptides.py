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
        ("data",  "d", "the path to the data folder containing the CSV tables")
    ]

    def initialize_options(self):
        self.force = False
        self.inplace = False
        self.data = None

    def finalize_options(self):
        _build_py = self.get_finalized_command("build_py")
        self.build_lib = _build_py.build_lib
        if self.data is None:
            self.data = os.path.relpath(os.path.join(__file__, "..", "peptides", "data"))

    def _load_tables(self):
        self.announce(f"loading data tables from {self.data!r}", level=2)
        tables = {}
        for group in os.listdir(self.data):
            if group == "__pycache__":
                continue
            if not os.path.isdir(os.path.join(self.data, group)):
                continue
            tables[group] = {}
            for filename in glob.glob(os.path.join(self.data, group, "*.csv")):
                member, _ = os.path.splitext(os.path.basename(filename))
                tables[group][member] = self._load_csv(filename)
        return tables

    def _load_csv(self, filename):
        self.announce(f"loading {filename!r}", level=1)
        with open(filename, "r") as f:
            return {row[0]:float(row[1]) for row in csv.reader(f)}

    def _generate_module(self, tables):
        n = sum(map(len, tables.values()))
        self.announce(f"building Python AST from {n!r} data tables", level=2)
        body = []
        for name, table in tables.items():
            stub = name.upper().replace(".", "_")
            assign_node = ast.Assign(
                targets=[ast.Name(id=stub, ctx=ast.Store())],
                value=ast.Dict(
                    keys=list(map(ast.Constant, table.keys())),
                    values=[
                        ast.Dict(
                            keys=list(map(ast.Constant, subtable.keys())),
                            values=list(map(ast.Constant, subtable.values())),
                        )
                        for subtable in table.values()
                    ]
                )
            )
            body.append(assign_node)
        return ast.Module(body=body)

    def _write_module(self, module, filename):
        self.announce(f"writing Python source to {filename!r}", level=2)
        with open(filename, "w") as f:
            f.write("# this file was automatically generated")
            f.write("# by the `python setup.py codegen` command\n")
            f.write("# DO NOT EDIT MANUALLY!\n")
            f.write(astor.to_source(module))
            f.write("\n")

    def run(self):
        # check astor is available
        if isinstance(astor, ImportError):
            raise RuntimeError("`astor` package is required for the `codegen` command") from astor
        # load data tables from R-formatted data file
        tables = self._load_tables()
        # generate a Python file containing the constants
        module_node = self._generate_module(tables)
        # write the data sources as Python
        output_file = os.path.join(self.build_lib, "peptides", "data", "tables.py")
        self.mkpath(os.path.dirname(output_file))
        self._write_module(module_node, output_file)
        # copy if inplace
        if self.inplace:
            library_file = os.path.join("peptides", "data", "tables.py")
            self.copy_file(output_file, library_file)


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
