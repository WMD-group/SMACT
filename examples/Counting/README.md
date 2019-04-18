There are 3 examples here:

- Generate_compositions_lists: A walkthrough of how to use `smact` to generate a list
of allowed compositions from a chosen search-space of elements. It shows you how to choose the elements 
you are interested in, and how to apply the standard smact_test with a certain stoichiometry threshold (see [docs here](https://smact.readthedocs.io/en/latest/examples.html#neutral-combinations) for more info).
It also shows you how to interface the output to [Pymatgen](http://pymatgen.org/) or [Pandas](https://pandas.pydata.org/). 

- Raw_combinations: Does not use `smact`, but uses itertools to calculate the raw number of possible 
element or species combinations for 2, 3 and 4-component compositions. 

- ElementCombinationsParallel.py: A script that reproduces the numbers in Table 1 of the 2016 Chem Paper 
[Computational Screening of All Stoichiometric Inorganic Materials](https://www.cell.com/chem/fulltext/S2451-9294(16)30155-3).
It is currently only Python2 compatible and is kept here as a record of the calculation methodology. 
