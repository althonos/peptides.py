#!/usr/bin/env python
# coding: utf-8

import ast
import csv
import glob
import math
import os

import setuptools
from distutils.command.build import build as _build

try:
    import astor
except ImportError as err:
    astor = err


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
        self._write_module(module_node, os.path.join("peptides", "data", "tables.py"))

if __name__ == "__main__":
    setuptools.setup(cmdclass=dict(codegen=codegen))
