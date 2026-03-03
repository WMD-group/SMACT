---
title: "SMACT v4: Enhanced Screening and Prediction Tools for Computational Materials Design"
tags:
  - Python
  - materials design
  - chemical heuristics
  - high-throughput screening
  - structure prediction
  - machine learning
authors:
  - name: Kinga O. Mastej
    orcid: 0009-0006-5656-6646
    affiliation: 1
  - name: Anthony Onwuli
    affiliation: 1
  - name: Alexander Moriarty
    affiliation: 2
  - name: Daniel Davies
    affiliation: 1
  - name: Ryan Nduma
    affiliation: 1
  - name: Keith Butler
    affiliation: 4
  - name: Jiwoo Lee
    affiliation: 3
  - name: Panyalak Detrattanawichai
    affiliation: 1
  - name: Alex Ganose
    affiliation: 5
  - name: Hyunsoo Park
    affiliation: 1
  - name: Rhys Goodall
    affiliation: 6
  - name: Tianshu Li
    affiliation: 1
  - name: Aron Walsh
    orcid: 0000-0001-5460-7033
    affiliation: 1
affiliations:
  - name: Department of Materials, Imperial College London, London, UK
    index: 1
  - name: Department of Chemical Engineering, University College London, London, UK
    index: 2
  - name: Department of Materials Science and Engineering, Yonsei University, Seoul, South Korea
    index: 3
  - name: Department of Chemistry, University College London, London, UK
    index: 4
  - name: Department of Chemistry, Imperial College London, London, UK
    index: 5
  - name: Radical AI, New York, USA
    index: 6
date: 3 March 2026
bibliography: paper.bib
---

# Summary

`SMACT` (Semiconducting Materials by Analogy and Chemical Theory) is a Python library for the high-throughput screening and discovery of inorganic materials. The name reflects the foundational analogy-based approaches to semiconductor prediction [@goodman1958; @pamplin1964]. Originally published as a JOSS paper in 2019 [@Davies2019], the package provided chemical filters based on charge neutrality, electronegativity, and elemental data to rapidly identify chemically plausible compositions from the vast space of possible element combinations [@davies_computational_2016]. Since then, SMACT has undergone substantial development spanning seven major releases (v2.1 through v4.0). The library now includes modules for crystal structure prediction via data-mined ionic substitution [@Hautier2011], dopant prediction, composition-to-property prediction using deep learning models [@Goodall2020], and tools for mapping and visualising inorganic crystal chemical space [@Park2025]. These additions transform SMACT from a compositional screening tool into a comprehensive platform for computational materials design.

# Statement of Need

The search space for new inorganic materials is immense---exceeding $4 \times 10^{12}$ quaternary element combinations [@davies_computational_2016]---and continues to grow as multi-component systems attract interest. Researchers need tools that can rapidly filter this space using established chemical principles before committing to expensive first-principles calculations or synthesis. The original SMACT addressed this need through charge neutrality and electronegativity screening. However, the field has since advanced: researchers now require not only compositional screening but also structure prediction, property estimation, and integration with modern machine learning workflows.

SMACT v4 meets these needs by providing a unified Python library where a researcher can (i) enumerate and filter compositions, (ii) predict likely crystal structures from data-mined substitution probabilities, (iii) estimate material properties directly from stoichiometry without requiring structural information, and (iv) identify candidate dopants for known host materials. The library is designed for desktop-scale computation, with data caching and multiprocessing support that enables screening of millions of compositions in minutes.

# State of the Field

Several established libraries serve the computational materials science community. `Pymatgen` [@ong2013] provides comprehensive tools for structure manipulation and thermodynamic analysis. The Atomic Simulation Environment (`ASE`) [@ase-paper] focuses on atomistic simulations. `Matminer` [@ward2018] specialises in featurisation and machine learning pipelines for materials data. These tools primarily operate on known structures or experimental data.

SMACT occupies a distinct niche: it works at the composition level, upstream of structure-based tools, to define which regions of chemical space merit further investigation. This top-down, composition-first approach complements structure-based workflows. The structure prediction and property prediction modules added since v2.0 provide natural bridges between SMACT's compositional screening and downstream structural analysis tools. The `ElementEmbeddings` integration [@Antunes2022] further connects SMACT to modern representation learning approaches for materials informatics [@Onwuli2024].

# Software Design

SMACT v4 retains the `Element` and `Species` classes at its core, providing access to tabulated data for all elements of the periodic table in arbitrary oxidation states and coordination environments. Building on this foundation, the following major modules have been added since the original publication.

**Structure prediction.** The `structure_prediction` subpackage implements data-mined ionic substitution methods [@Hautier2011]. A `CationMutator` class computes substitution probabilities from lambda tables derived from a statistical analysis of known compounds. The `StructurePredictor` applies these probabilities to predict likely crystal structures for new compositions by analogy with parent structures stored in a local SQLite database (`StructureDB`), which can be populated from the Materials Project.

**Dopant prediction.** The `dopant_prediction` module (`Doper` class) enables high-throughput identification of p-type and n-type dopants for a given host composition. It combines oxidation state filters with ionic radius constraints and optional species embedding similarities to rank candidate dopants by substitution probability and site selectivity.

**Property prediction.** The `property_prediction` subpackage provides composition-to-property prediction using pre-trained ROOST (Representation Learning from Stoichiometry) models [@Goodall2020]. Users can predict band gaps directly from chemical formulae, with the extensible architecture designed to support additional properties as pre-trained models become available. Ensemble-based uncertainty quantification is included, and a model registry and caching system manage pre-trained checkpoints.

**Oxidation states and metallicity.** The `oxidation_states` module implements a probabilistic model for predicting the likelihood of metal species coexisting in compounds based on anion-dependent statistics [@davies2018]. SMACT v3 updated the default oxidation state data from the original SMACT14 set to ICSD24, reflecting more recent experimental observations. The `metallicity` module provides scoring functions based on electronegativity differences to distinguish ionic compounds from intermetallic phases.

**Crystal space utilities.** The `crystal_space` subpackage provides tools for systematically mapping inorganic composition space [@Park2025], including Materials Project data acquisition, SMACT-based composition generation, and dimensionality-reduced visualisation of compositional embeddings.

**Screening enhancements.** The core `smact_filter` and `smact_validity` functions have been extended with support for mixed-valence compositions, switchable oxidation state datasets (ICSD24, ICSD16, SMACT14, Pymatgen, Wikipedia), early short-circuiting for performance, and consensus-based oxidation state filtering.

**Code quality.** SMACT v4 enforces full type annotations checked by `pyright` in standard mode, linting with `ruff` across all rule categories, and maintains greater than 85% test coverage across 318 tests running on three Python versions and three operating systems via GitHub Actions.

# Research Impact

SMACT has been used in a number of published research studies. Park et al. [@Park2025] used SMACT's screening and crystal space tools to map the landscape of inorganic crystal chemistry. Onwuli et al. [@Onwuli2024] employed SMACT's species representations for materials informatics, developing ionic embeddings that capture oxidation-state-dependent chemical behaviour. The structure prediction module has supported studies on data-driven discovery of new compounds using ionic substitution methods. The dopant prediction capabilities have been applied to identify candidate dopants for semiconductor and energy materials. Park et al. [@Park2025synthesis] used SMACT in a workflow for closing the synthesis gap in computational materials design. Nduma et al. [@Nduma2025crystalyse] integrated SMACT into the Crystalyse multi-tool agent for autonomous materials design.

# AI Usage Disclosure

Generative AI tools (Claude, Anthropic) were used to assist with code modernisation tasks during the v4.0.0 development cycle, including type annotation additions, test migration, and documentation formatting. All AI-generated code was reviewed and validated by the authors. The first draft of this manuscript was prepared with AI assistance and reviewed by all co-authors.

# Author Contributions

[KOM](https://github.com/KingaMas) is the current maintainer. She contributed consensus-based oxidation state filtering for `smact_validity` and `smact_filter`, and led the v4.0.0 release including mixed-valence support, a deep code audit, full type annotations, security hardening, test migration to pytest with coverage raised to over 85%, and CI consolidation.
[AO](https://github.com/AntObi) was the primary developer from v2.2 to v3.1, contributing the dopant prediction module, oxidation states improvements, crystal space utilities, and extensive maintenance.
[AM](https://github.com/a-ws-m) designed and implemented the structure prediction subpackage and benchmarking framework.
[DWD](https://github.com/dandavies99) is the original lead developer and contributed to CI infrastructure and code quality tooling.
[RN](https://github.com/ryannduma) implemented the property prediction subpackage, metallicity module, and screening performance optimisations.
[KTB](https://github.com/keeeto) is an original author who continued contributing to the codebase.
[JL](https://github.com/JiwooChloeLee) contributed the initial dopant prediction implementation and charge state comparison fixes.
[PD](https://github.com/Panyalak) contributed to lattice parameter calculations.
[AMG](https://github.com/utf) contributed to code quality and structure.
[HP](https://github.com/hspark1212) contributed the crystal space visualisation tools and element data corrections.
[REAG](https://github.com/comprhys) fixed a screening bug related to electronegativity handling.
[TL](https://github.com/lits19) contributed to code reorganisation and utility refactoring.
[AW](https://wmd-group.github.io/) supervised the project and contributed to the codebase.

# Acknowledgements

We acknowledge contributions from all members of the Walsh Materials Design group who have provided feedback, testing, and feature requests. We thank the JOSS editors and reviewers of the original 2019 paper. This work was supported by the Engineering and Physical Sciences Research Council (EPSRC), UK.

# References
