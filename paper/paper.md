---
title: 'mol-ellipsize: A Python package for molecular size calculation'
tags:
  - Python
  - rdkit
  - chemistry
  - cheminformatics
  - molecules
authors:
  - name: Andrew Tarzia*
    orcid: 0000-0001-8797-8666
    affiliation: "1, 2"
  - name: Christian Doonan
    orcid: 0000-0003-2822-0956
    affiliation: 2
  - name: David Huang*
    orcid: 0000-0003-2048-4500
    affiliation: 2
affiliations:
 - name: Department of Chemistry, Molecular Sciences Research Hub, Imperial College London, White City Campus, Wood Lane, London, W12 0BZ, UK
   index: 1
 - name: Department of Chemistry and Centre for Advanced Nanomaterials, The University of Adelaide, South Australia 5005, Australia
   index: 2
date: 18 June 2021
bibliography: paper.bib

---

# Summary

`mol-ellipsize` is a Python package for calculating the size of a molecule based on the fit of an ellipsoid to its shape. The shape can be arbitrarily defined, but a default, based on the molecule's van der Waals (vdW) volume, is implemented and agrees well with reported kinetic diameters for a wide variety of small molecules. `mol-ellipsize` provides a workflow for automating the calculation of molecular size using common cheminformatics tools with an ASCII-string representation of the molecule's chemical structure as the only input, such that minimal human intervention is required. By employing an ellipsoid fit to the molecular volume of automatically generated conformers, `mol-ellipsize` enables more straightforward and accurate consideration (compared with methods based on spherical or cylindrical measures) of the effects of molecular shape and conformational flexibility on estimates of molecular size relevant to mass transport in porous materials.

# Statement of need

The size and shape of a molecule determines many of its properties.
Specifically, in the design of molecular adsorbents or sieves, the size of the target molecule (or guest) is crucial. Materials chemists often use the size of guest molecules as descriptors for screening potential materials for whether certain guests will diffuse through them. As molecules increase in complexity, the calculation of their size becomes difficult. Additionally, the manual calculation of the size of many different molecules or molecular conformers is intractable. `mol-ellipsize` solves both of these issues by providing an efficient, automated workflow for calculating the size of a molecule and its conformers. After describing the algorithm, we demonstrate a research application of `mol-ellipsize` to mass transport in porous functional materials.

# Calculation of molecular size

In `mol-ellipsize`, we have implemented an algorithm to calculate the size of a molecule from its chemical graph. The algorithm effectively calculates the cross-sectional size of the van der Waals (vdW) volume of a given molecule in a similar fashion to previous approaches used to determine the critical diameter of an adsorbate (the critical diameter is derived from a cylindrical fit to the vdW volume of a molecule [@webster1998; @polyukhov2019]). However, `mol-ellipsize` automates this process and uses an ellipsoid to fit the vdW cloud rather than a cylinder. Importantly, `mol-ellipsize` provides a low-cost automated conformer generation algorithm, such that flexibility can be considered in a user-friendly way.

\autoref{fig:size_method} shows an example workflow of `mol-ellipsize` usage. We use RDKit (a cheminformatics Python toolkit [@landrum]) to generate conformers [@riniker2015] of a molecule and to calculate the vdW volume of each of those conformers, which are input to the ellipsoid fitting algorithm. The vdW volume definition is based on a grid, which is defined by the ``box margin`` and ``grid spacing`` input arguments that specify the extent of the grid beyond the edges of the molecule in each dimension and the resolution of the grid, respectively. `mol-ellipsize` calculates the minimum-volume enclosing ellipsoid (shown schematically in black in \autoref{fig:size_method}b) for the points in the vdW volume of each conformer using a minimization algorithm based on the Khachiyan algorithm (we used a tolerance of 0.1 for the relative change in the solution as a stopping criterion).[@moshtagh; @imelfort] An ellipsoid is defined by its three principal diameters, with its second largest (or intermediate) diameter defining the smallest sized cross-section that is required to diffuse through the pores of a porous material (assuming cylindrical pores [@webster1998]). It is then possible to assign the size of a molecule based on the minimum intermediate diameter ($d$) of all of its conformers. We have parameterised this algorithm for a series of molecules (\autoref{tbl:1}), such that default settings in `mol-ellipsize` achieves molecular sizes that agree with kinetic diameters. The kinetic diameter is a length scale characterising the intermolecular separation of gas-phase collisions [@zeolite-molecular-sieves] and is determined experimentally from second virial coefficient or gas viscosity data [@matteucci2006]. It is strongly correlated with the diffusion coefficient in porous media and is thus commonly used to predict rates of mass transport in porous materials [@zhang2012; @zhang2013]. \autoref{fig:size_method_parity}a shows good agreement between the calculated molecular size $d$ and reported kinetic diameters for all molecules in \autoref{tbl:1} using a box margin of 4 $\overset{\circ}{\mathrm{A}}$, a grid spacing of 0.5 $\overset{\circ}{\mathrm{A}}$, $N=100$, and a vdW scale parameter (the value by which the vdW parameters are scaled) of 0.9. We find that the results are not very sensitive to parameter choice for a physically reasonable parameter range. Users should consider the trade-off between accuracy and computational efficiency. However, the methodology appears sufficient to approximate the kinetic diameters of small molecules accurately.


![(a) Sequence of steps in an example calculation of the molecular size of $n$-octane from its SMILES string (an ASCII representation of the molecule). Multiple 3D conformers (100 in this work; a subset is shown in distinct colors in (b) along with the distribution of all diameters) are generated and the minimum-volume enclosing ellipsoid (shown schematically as black dashed lines in (b)) that encompasses a grid representation of the vdW volume (colored surfaces in (b)) of each conformer is calculated. The molecular size $d$ of a molecule is given by the minimum intermediate diameter (blue line in (c)) of all of its conformer's ellipsoids.\label{fig:size_method}](size_method.pdf)


![Parity plots of the calculated molecular size $d$ versus the (a) reported kinetic diameters for all molecules in \autoref{tbl:1} and (b) critical diameters of small molecules extracted from ref. @webster1998.\label{fig:size_method_parity}](main_parities.pdf)

|    name          | kinetic diameter |          name          |               kinetic diameter              |
|:----------------:|:----------------:|:----------------------:|:-------------------------------------------:|
|     He           |       2.551      |     dimethyl ether     |                    4.307                    |
|     Ne           |       2.82       |         ethane         |                    4.443                    |
|     Ar           |       3.542      |         ethene         |                    4.163                    |
|     Kr           |       3.655      |         ethanol        |                    4.530                    |
|     Xe           |       4.047      |       $n$-propane      |                  4.3--5.118                 |
| H$_2$            | 2.827--2.89      | cyclopropane           | 4.23--4.807                                 |
| Cl$_2$           | 4.217            | propene                | 4.678                                       |
| Br$_2$           | 4.296            | acetone                | 4.600                                       |
| CO$_2$           | 3.3              | $n$-butane             | 4.687                                       |
| O$_2$            | 3.467            | 1-butene               | 4.5                                         |
| N$_2$            | 3.64--3.80       | $i$-butane             | 5.278                                       |
| H$_2$O           | 2.641            | 2,2-dimethylbutane     | 6.2                                         |
| NO               | 3.492            | cis-2-butene           | 4.23                                        |
| CO               | 3.69             | 1,3-butadiene          | 5.2                                         |
| N$_2$O           | 3.828            | $n$-pentane            | 4.5                                         |
| HCl              | 3.339            | $i$-pentane            | 5                                           |
| HBr              | 3.353            | neo-pentane            | 6.2--6.464                                  |
| CS$_2$           | 4.483            | 2-methyl pentane       | 5.5                                         |
| COS              | 4.130            | 2,2,4-trimethylpentane | 6.2                                         |
| SO$_2$           | 4.112            | 3-methylpentane        | 5.5                                         |
| H$_2$S           | 3.623            | $n$-hexane             | 4.3                                         |
| NH$_3$           | 2.900            | $n$-heptane            | 4.3                                         |
| NF$_3$           | 3.62             | $n$-octane             | 4.3                                         |
| CCl$_2$F$_2$     | 5.0              | cyclohexane            | 6--6.182                                    |
| CH$_3$Cl         | 4.182            | benzene                | 5.349--5.85                                 |
| CH$_2$Cl$_2$     | 4.898            | ethyl-benzene          | 5.8                                         |
| CHCl$_3$         | 5.389            | $p$-xylene             | 5.8                                         |
| CCl$_4$          | 5.947            | $m$-xylene             | 6.8                                         |
| CF$_4$           | 4.662            | $o$-xylene             | 6.8                                         |
| C$_2$F$_6$       | 5.1              | $i$-butene             | 4.8[@zhang2012]                             |
| $n$-C$_6$F$_{14}$| 7                | 1-butanol              | 4.5[@zhang2013]                             |
| methane          | 3.758            | 2,3-dimethylbutane     | 5.6[@zhang2013]                             |
| methanol         | 3.626            | 1,2,4-trimethylbenzene | 7.6[@zhang2013]                             |
| acetylene        | 3.3              | mesitylene             | 8.2\textsuperscript{\emph{a}}[@webster1998] |
| toluene          | 5.25             |                        |                                             |

Table:  Kinetic diameters of molecules (in $\overset{\circ}{\mathrm{A}}$) used to parameterize our methodology for calculating the molecular size. All kinetic diameters were taken from ref. @li2009 unless otherwise cited. Where applicable, the smaller value of a range was used in \autoref{fig:size_method_parity}\label{tbl:1}.


# Application to research

Enzymes are a class of protein whose function is to catalyze the biochemical transformation of a substrate to a product. Improved durability is a crucial step toward the commercial application of enzymes as industrial catalysts.[@schmid2001] Through encapsulation in metal--organic frameworks (MOFs) matrices, enzymes are able to retain their activity in harsh conditions (e.g. elevated temperatures or proteolytic media) due to the protection afforded by the MOF matrix.[@liang2015] Additionally, the MOF matrix is expected to afford size-selective transport of substrates to the active site of an enzyme via its pore network. The cost of purified enzymes and the time required to prepare and screen new enzymes is significant. Hence, accurate predictions of which reactions may be possible inside a given material before doing costly experiments are valuable. Therefore, we implemented `mol-ellipsize` to enable efficient screening of large numbers of reactants and products of enzymatic reactions for their ability to diffuse within MOFs, based on their molecular size. We used the KEGG database[@kanehisa2000; @kanehisa2016; @kanehisa2017] of ~12000 curated enzymatic reactions to compile a list of possible reactions with known molecular components. From ~12000 KEGG reactions we extracted a chemical space of 5640 unique molecules (that were readable by RDKit with a mass limit of 500 g/mol; a mass limit was applied because larger molecules are not expected to diffuse through typical MOFs). \autoref{fig:chem_space} shows the calculated minimum intermediate diameter, $d$, of these molecules as a function of the number of heavy atoms in the molecule. Unsurprisingly, there is a correlation between size and the number of heavy atoms, which begins to plateau as larger molecules tend to be more flexible and are, therefore, able to adopt conformations with smaller intermediate diameters than would be predicted from a naive consideration of their chemical formulae. \autoref{fig:chem_space} clearly shows that while the majority of the chemical space is too large to diffuse through ZIF-8 (a prototypical porous material), there are still a significant number of molecules (~1000) with $d<$ 6.6 $\overset{\circ}{\mathrm{A}}$ (the approximate diffusion threshold for ZIF-8[@zhang2012; @zhang2013; @verploegh2015; @diestel2012; @zhu2016]), which should be able to diffuse through ZIF-8. Therefore, we would expect that novel candidate enzyme reactions can be found. Furthermore, the generalisability of `mol-ellipsize` allows this approach to expand that chemical space through synthetic modifications to the framework or reaction components.


![Distribution of the minimum intermediate diameter $d$ calculated by `mol-ellipsize} as a function of the number of heavy atoms in each molecule of all collected molecules. The shaded region in indicates the approximate range for the threshold for diffusion through ZIF-8 from the literature.\label{fig:chem_space}](joss_chemspace.pdf)

All code used to extract molecules from KEGG is available at [github.com/andrewtarzia/enzyme_screen](https://github.com/andrewtarzia/enzyme_screen.git).


# Acknowledgements

A. T. was supported by an Australian Government RTP Scholarship and a PhD top-up scholarship from CSIRO Division of Materials Science and Engineering. A. T. would like to thank Dr Jesse Teo, Dr Natasha Maddigan, Dr Weibin Liang, Dr Justin Spence, A/Prof Stephen Bell, Dr Austin Mroz and Dr Kim E. Jelfs for discussions about this work.

# References

