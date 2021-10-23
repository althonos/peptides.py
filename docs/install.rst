Installation
============

PyPi
^^^^

``peptides.py`` is hosted on GitHub, but the easiest way to install it is to
download the latest release from its `PyPi repository <https://pypi.python.org/pypi/pyrodigal>`_.

.. code:: console

	$ pip install --user pyrodigal

.. Conda
.. ^^^^^
..
.. Pronto is also available as a `recipe <https://anaconda.org/bioconda/pyrodigal>`_
.. in the `bioconda <https://bioconda.github.io/>`_ channel. To install, simply
.. use the ``conda`` installer:
..
.. .. code:: console
..
.. 	 $ conda install -c bioconda pyrodigal


GitHub + ``pip``
^^^^^^^^^^^^^^^^

If, for any reason, you prefer to download the library from GitHub, you can clone
the repository and install the repository by running (with the admin rights):

.. code:: console

  $ git clone --recursive https://github.com/althonos/peptides.py
	$ pip install --user ./peptides.py

.. caution::

    Keep in mind this will install always try to install the latest commit,
    which may not even build, so consider using a versioned release instead.


GitHub + ``setuptools``
^^^^^^^^^^^^^^^^^^^^^^^

If you do not want to use ``pip``, you can still clone the repository and
run the ``setup.py`` file manually, although you will need to install the
build dependencies (mainly `astor <https://pypi.org/project/astor>`_):

.. code:: console

	$ git clone --recursive https://github.com/althonos/peptides.py
	$ cd peptides.py
	# python setup.py build_py install

.. Danger::

    Installing packages without ``pip`` is strongly discouraged, as they can
    only be uninstalled manually, and may damage your system.
