Installation
------------

Linkapy requires Python >= 3.10, < 3.13.

Linkapy can be installed via pypi:

.. code-block:: bash

    pip install linkapy

A development version can be installed from github, by cloning the repository.
Note that a development version requires you to have `Maturin <https://github.com/PyO3/maturin>`_ installed.
Additionaly, you need to have `Rust <https://www.rust-lang.org/tools/install>`_ (1.85 or newer, for 2024 edition support) installed as well.
Linkapy can then be installed with:

.. code-block:: bash

    git clone git@github.com:WardDeb/linkapy.git
    cd linkapy
    pip install .


Can also be done with maturin:

.. code-block:: bash

    git clone git@github.com:WardDeb/linkapy.git
    cd linkapy
    maturin develop --release

There are additional dependencies for documentation, building, and developing/testing. These can be included too:

.. code-block:: bash

    maturin develop --release --extras docs,dev,build

or

.. code-block:: bash

    pip install .[docs,dev,build]

Linkapy can also be used through pixi. After cloning the repository, you can install / use linkapy with:

.. code-block:: bash

    pixi run linkapy -h

