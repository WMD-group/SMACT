[![DOI](https://joss.theoj.org/papers/10.21105/joss.01361/status.svg)](https://doi.org/10.21105/joss.01361)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5553202.svg)](https://doi.org/10.5281/zenodo.5553202)
[![Documentation Status](https://readthedocs.org/projects/smact/badge/?version=latest)](http://smact.readthedocs.org/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![python version](https://img.shields.io/pypi/pyversions/smact)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPi](https://img.shields.io/pypi/v/smact)](https://pypi.org/project/SMACT/)
[![Conda](https://anaconda.org/conda-forge/smact/badges/version.svg)](https://anaconda.org/conda-forge/smact)
[![GitHub issues](https://img.shields.io/github/issues-raw/WMD-Group/SMACT)](https://github.com/WMD-group/SMACT/issues)
![dependencies](https://img.shields.io/librariesio/release/pypi/smact)
[![CI Status](https://github.com/WMD-group/SMACT/actions/workflows/ci.yml/badge.svg)](https://github.com/WMD-group/SMACT/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/WMD-group/SMACT/branch/master/graph/badge.svg?token=UtgVxjoYNP)](https://codecov.io/gh/WMD-group/SMACT)
![PyPI - Downloads](https://img.shields.io/pypi/dm/smact)

# SMACT

**Semiconducting Materials from Analogy and Chemical Theory** (SMACT) is a collection of rapid screening and informatics tools that uses data about chemical elements.

- **Documentation:** <https://smact.readthedocs.io/en/latest/>
- **Examples:** <https://smact.readthedocs.io/en/latest/examples.html>

![A blue interface with the text "SMACT v3" at the top. Below that, there is a label "Materials Search" followed by two radio buttons: "Hi-fi" and "Lo-fi". The "Lo-fi" button is currently selected](SMACT.png)

_If you torture the data enough, nature will always confess_ - Roland Coase (from 'How should economists choose?')

## Statement of need

There is a strong demand for functional materials across a wide range of technologies. The motivation can include cost reduction, performance enhancement, or to enable a new application. We have developed low-cost procedures for screening hypothetical materials. This framework can be used for simple calculations on your own computer. SMACT follows a top-down approach where a set of element combinations is generated and then screened using rapid chemical filters. It can be used as part of a multi-technique workflow or to feed artificial intelligence models for materials.

![A gif depicting using the SMACT code. The first lines of code show how SMACT can be used to access properties of Iron (Fe) by create an Fe Element object and then accessing the oxidation states, pauling electronegativity. The next line after these shows the use of the smact_filter function for the Fe-Cu-O chemical system followed by the lists of possible compositions. ](smact_simple.gif)

## Getting started

Features are accessed through Python scripts, importing classes and functions as needed.
The best place to start is looking at [the docs](https://smact.readthedocs.io/en/latest/), which highlight some simple examples of how these classes and functions can be usede
Use cases are available in our [examples](https://smact.readthedocs.io/en/latest/examples.html) and [tutorials](https://smact.readthedocs.io/en/latest/tutorials.html) folders.

## Code features

- At the core of SMACT are [Element](https://smact.readthedocs.io/en/latest/smact.html#smact.Element) and [Species](https://smact.readthedocs.io/en/latest/smact.html#smact.Species) (element in a given oxidation state) classes that have various properties associated with them.

- Oxidation states that are accessible to each element are included in their properties.

- Element compositions can be screened through based on the heuristic filters of charge neutrality and electronegativity order. This is handled using the [screening module](https://smact.readthedocs.io/en/latest/smact.screening.html) and [this publication](<https://www.cell.com/chem/fulltext/S2451-9294(16)30155-3>) describes the underlying theory. An example procedure is [outlined in the docs](https://smact.readthedocs.io/en/latest/examples/filter.html).

- Further filters can be applied to generated lists of compositions in order to screen for particular properties. These properties are either intrinsic properties of elements or are calculated for compositions using the [properties module](https://smact.readthedocs.io/en/latest/smact.properties.html). For example:

  - An application is shown in [this publication](https://pubs.rsc.org/en/content/articlehtml/2018/sc/c7sc03961a), in which 160,000 chemical compositions are screened based on optical band gap calculated using the [solid-state energy scale](https://www.sciencedirect.com/science/article/pii/S0022459615300888).
  - The [oxidation_states module](https://smact.readthedocs.io/en/latest/smact.oxidation_states.html) can be used to filter out compositions containing metals in unlikely oxidation states according to [a data-driven model](https://pubs.rsc.org/en/content/articlelanding/2018/fd/c8fd00032h#!divAbstract).

- Compositions can also be filtered based on sustainability via the abundance of elements in the Earth's crust or via the [HHI scale](https://pubs.acs.org/doi/10.1021/cm400893e).

- Compositions can be converted for use in Pymatgen or for representation to machine learning algorithms ([see this example](https://smact.readthedocs.io/en/latest/tutorials/smact_generation_of_solar_oxides.html)) and the related [ElementEmbeddings](https://github.com/WMD-group/ElementEmbeddings) package.

- The code also has tools for manipulating common crystal lattice types:
  - Certain structure types can be built using the [builder module](https://smact.readthedocs.io/en/latest/smact.builder.html)
  - Lattice parameters can be estimated using ionic radii of the elements for various common crystal structure types using the [lattice_parameters module](https://smact.readthedocs.io/en/latest/smact.lattice_parameters.html).
  - The [lattice module](https://smact.readthedocs.io/en/latest/smact.lattice.html) and [distorter module](https://smact.readthedocs.io/en/latest/smact.distorter.html) rely on the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/) and can be used to generate unique atomic substitutions on a given crystal structure.
  - The [structure prediction](https://smact.readthedocs.io/en/latest/smact.structure_prediction.html) module can be used to predict the structure of hypothetical compositions using species similarity measures.
  - The [dopant prediction](https://smact.readthedocs.io/en/latest/smact.dopant_prediction.html) module can be used to facilitate high-throughput predictions of p-type and n-type dopants of multicomponent solids.

## List of modules

- **smact** library containing:
  - **\_\_init\_\_.py** Contains the core `Element` and `Species` classes.
  - **data_loader.py** Handles the loading of external data used to initialise the core `smact.Element` and `smact.Species` classes.
  - **screening.py** Used for generating and applying filters to compositional search spaces.
  - **properties.py** A collection of tools for estimating useful properties based on composition.
  - **lattice.py** Given the sites, multiplicities and possible oxidation states
    at those sites, this reads from the database and generates all possible
    stoichiometries.
  - **builder.py** Builds some common lattice structures, given the chemical
    composition.
  - **lattice_parameters.py** Estimation of lattice parameters for various lattice types using covalent/ionic radii.
  - **distorter.py** A collection of functions for enumerating and then
    substituting on inequivalent sites of a sub-lattice.
  - **oxidation_states.py**: Used for predicting the likelihood of species coexisting in a compound based on a statistical model.
  - **structure_prediction**: A submodule which contains a collection of tools for facilitating crystal structure predictions via ionic substitutions
  - **dopant_prediction**: A submodule which contains a collections of tools for predicting dopants.
  - **utils.py** A collection of utility functions used throughout the codebase.

## Requirements

The main language is Python 3 and has been tested using Python 3.10+.
Basic requirements are Numpy and Scipy.
The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase) (ASE), [spglib](http://atztogo.github.io/spglib), and [pymatgen](https://pymatgen.org) are also required for many components.

## Installation

The latest stable release can be installed via pip which will automatically set up other Python packages as required:

    pip install smact

Optional dependencies can also be installed. These enable full replication of the examples and tutorials

    pip install "smact[optional]"

SMACT is also available via conda through the conda-forge channel on Anaconda Cloud:

    conda install -c conda-forge smact

Alternatively, the very latest version can be installed using:

    pip install git+https://github.com/WMD-group/SMACT.git

For developer installation SMACT can be installed from a copy of the source
repository (<https://github.com/wmd-group/smact>); this will be preferred if using experimental code branches.

To clone the project from GitHub and make a local installation:

    git clone https://github.com/wmd-group/smact.git
    cd smact
    pip install --user -e .

With -e pip will create links to the source folder so that that changes
to the code will be immediately reflected on the PATH.

## License and attribution

Python code and original data tables are licensed under the MIT License.

## Development notes

### Bugs, features and questions

Please use the [Issue Tracker](https://github.com/WMD-group/smact/issues) to report bugs or request features in the first instance. While we hope that most questions can be answered by searching [the docs](https://smact.readthedocs.io/en/latest/), we welcome new questions on the issue tracker, especially if they helps us improve the docs! For other queries about any aspect of the code, please contact Anthony Onwuli (maintainer) by [e-mail](mailto:anthony.onwuli16@imperial.ac.uk).

### Code contributions

We are always looking for ways to make SMACT better and more useful to the wider community; contributions are welcome. Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

- Code style should comply with [PEP8](http://www.python.org/dev/peps/pep-0008) where possible. [Google's house style](https://google.github.io/styleguide/pyguide.html) is also helpful, including a good model for docstrings.
- Please use comments liberally when adding nontrivial features, and take the chance to clean up other people's code while looking at it.
- Add tests wherever possible, and use the test suite to check if you broke anything.
- Look at the [contributing guide](CONTRIBUTING.md) for more information.

### Tests

We use integrated testing on GitHub via [GitHub Actions](https://github.com/features/actions). Testing modules should be pass/fail and wrapped into **tests/test_core.py** or another **tests/test_something.py** file added, if appropriate.
Run the tests using `python -m pytest -v`.(The final `-v` is optional and adds more detail to the output.)

## References

[H. Park et al.,
"Mapping inorganic crystal chemical space" _Faraday Discuss._ (2024)](https://pubs.rsc.org/en/content/articlelanding/2024/fd/d4fd00063c)

[D. W. Davies et al.,
"SMACT: Semiconducting Materials by Analogy and Chemical Theory" _JOSS_ **4**, 1361 (2019)](https://joss.theoj.org/papers/7efd2f2ad60d25bdccee3fbd3fc11448)

[D. W. Davies et al.,
"Materials discovery by chemical analogy: role of oxidation states in structure prediction" _Faraday Discuss._ **211**, 553 (2018)](https://pubs.rsc.org/en/Content/ArticleLanding/2018/FD/C8FD00032H)

[D. W. Davies et al.,
"Computational screening of all stoichiometric inorganic materials" _Chem_ **1**, 617 (2016)](<http://www.cell.com/chem/abstract/S2451-9294(16)30155-3>)

[B. R. Pamplin, "A systematic method of deriving new semiconducting compounds by structural analogy", _J. Phys. Chem. Solids_
**25**, 675 (1964)](http://www.sciencedirect.com/science/article/pii/0022369764901763)
