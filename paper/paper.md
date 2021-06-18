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

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }


# Calculation of molecular size

In `mol-ellipsize`, we have implemented an algorithm to calculate the size of a molecule from its chemical graph. The algorithm effectively calculates the cross-sectional size of the van der Waals (vdW) volume of a given molecule in a similar fashion to previous approaches used to determine the critical diameter of an adsorbate (the critical diameter is derived from a cylindrical fit to the vdW volume of a molecule)[@webster1998,@polyukhov2019] However, `mol-ellipsize` automates this process and uses an ellipsoid to fit the vdW cloud rather than a cylinder. Importantly, `mol-ellipsize` provides a low-cost automated conformer generation algorithm, such that flexibility can be considered in a user-friendly way.

Figure~\ref{fgr:size_method} shows an example workflow of `mol-ellipsize` usage. We use RDKit (a cheminformatics Python toolkit)[@landrum] to generate conformers[@riniker2015] of a molecule and to calculate the vdW volume of each of those conformers, which are input to the ellipsoid fitting algorithm. The vdW volume definition is based on a grid, which is defined by the `box margin' and `grid spacing' input arguments that specify the extent of the grid beyond the edges of the molecule in each dimension and the resolution of the grid, respectively. `mol-ellipsize` calculates the minimum-volume enclosing ellipsoid (shown schematically in black in Figure~\ref{fgr:size_method}b) for the points in the vdW volume of each conformer using a minimization algorithm based on the Khachiyan algorithm (we used a tolerance of 0.1 for the relative change in the solution as a stopping criterion).[@moshtagh; @imelfort] An ellipsoid is defined by its three principal diameters, with its second largest (or intermediate) diameter defining the smallest sized cross-section that is required to diffuse through the pores of a porous material (assuming cylindrical pores).[@webster1998] It is then possible to assign the size of a molecule based on the minimum intermediate diameter (`d`) of all of its conformers. We have parameterised this algorithm for a series of molecules (Table~\ref{stbl:test_molecules1}), such that default settings in `mol-ellipsize` achieves molecular sizes that agree with kinetic diameters. The kinetic diameter is a length scale characterising the intermolecular separation of gas-phase collisions[@zeolite-molecular-sieves] and is determined experimentally from second virial coefficient or gas viscosity data.[@matteucci2006] It is strongly correlated with the diffusion coefficient in porous media and is thus commonly used to predict rates of mass transport in porous materials.[@zhang2012; @zhang2013] Figure~\ref{fgr:size_method_parity}a shows good agreement between the calculated molecular size $d$ and reported kinetic diameters for all molecules in Table~\ref{stbl:test_molecules1} using a box margin of \SI{4}{\angstrom}, a grid spacing of \SI{0.5}{\angstrom}, $N=100$, and a vdW scale parameter (the value by which the vdW parameters are scaled) of 0.9. We find that the results are not very sensitive to parameter choice for a physically reasonable parameter range. Users should consider the trade-off between accuracy and computational efficiency. However, the methodology appears sufficient to approximate the kinetic diameters of small molecules accurately.

\begin{figure}[ht!]
\centering
\includegraphics{size_method}
\caption{
(a) Sequence of steps in an example calculation of the molecular size of `n}-octane from its SMILES string (an ASCII representation of the molecule).
Multiple 3D conformers (100 in this work; a subset is shown in distinct colors in (b) along with the distribution of all diameters) are generated and the minimum-volume enclosing ellipsoid (shown schematically as black dashed lines in (b)) that encompasses a grid representation of the vdW volume (colored surfaces in (b)) of each conformer is calculated.
The molecular size $d$ of a molecule is given by the minimum intermediate diameter (blue line in (c)) of all of its conformer's ellipsoids.
}
\label{fgr:size_method}
\end{figure}

\begin{figure}[t!]
    \centering
    \includegraphics{main_parities}
    \caption{
        Parity plots of the calculated molecular size $d$ versus the (a) reported kinetic diameters for all molecules in Table~\ref{stbl:test_molecules1} and (b) critical diameters of small molecules extracted from ref. @webster1998.
    }
    \label{fgr:size_method_parity}
\end{figure}


\begin{table}[hp!]
    \centering
    \caption{
        Kinetic diameters of molecules used to parameterize our methodology for calculating the molecular size. All kinetic diameters were taken from ref. @li2009 unless otherwise cited. Where applicable, the smaller value of a range was used in Figure~\ref{fgr:size_method_parity}.
    }
    \begin{threeparttable}
        \begin{tabular}{lclc}
            \hline
            name & kinetic diameter [\AA] &name & kinetic diameter [\AA] \\
            \hline
            \ce{He}&2.551&dimethyl ether&4.307 \\
            \ce{Ne}&2.82&ethane&4.443 \\
            \ce{Ar}&3.542&ethene&4.163 \\
            \ce{Kr}&3.655&ethanol&4.530 \\
            \ce{Xe}&4.047&`n}-propane&4.3--5.118 \\
            \ce{H2}&2.827--2.89&cyclopropane&4.23--4.807 \\
            \ce{Cl2}&4.217&propene&4.678 \\
            \ce{Br2}&4.296&acetone&4.600 \\
            \ce{CO2}&3.3&`n}-butane&4.687 \\
            \ce{O2}&3.467&1-butene&4.5 \\
            \ce{N2}&3.64--3.80&`i}-butane&5.278 \\
            \ce{H2O}&2.641&2,2-dimethylbutane&6.2 \\
            \ce{NO}&3.492&cis-2-butene&4.23 \\
            \ce{CO}&3.69&1,3-butadiene&5.2 \\
            \ce{N2O}&3.828&`n}-pentane&4.5 \\
            \ce{HCl}&3.339&`i}-pentane&5 \\
            \ce{HBr}&3.353&neo-pentane&6.2--6.464 \\
            \ce{CS2}&4.483&2-methyl pentane&5.5 \\
            \ce{COS}&4.130&2,2,4-trimethylpentane&6.2 \\
            \ce{SO2}&4.112&3-methylpentane&5.5 \\
            \ce{H2S}&3.623&`n}-hexane&4.3 \\
            \ce{NH3}&2.900&`n}-heptane&4.3 \\
            \ce{NF3}&3.62&`n}-octane&4.3 \\
            %			\ce{SF6}&5.128&toluene&5.25 \\
            \ce{CCl2F2}&5.0&cyclohexane&6--6.182 \\
            \ce{CH3Cl}&4.182&benzene&5.349--5.85 \\
            \ce{CH2Cl2}&4.898&ethyl-benzene&5.8 \\
            \ce{CHCl3}&5.389&`p}-xylene&5.8 \\
            \ce{CCl4}&5.947&`m}-xylene&6.8 \\
            \ce{CF4}&4.662&`o}-xylene&6.8 \\
            \ce{C2F6}&5.1&`i}-butene&4.8[@zhang2012] \\
            `n}-\ce{C6F14}&7&1-butanol&4.5[@zhang2013] \\
            methane&3.758&2,3-dimethylbutane&5.6[@zhang2013] \\
            methanol&3.626&1,2,4-trimethylbenzene&7.6[@zhang2013] \\
            acetylene&3.3&mesitylene&8.2\textsuperscript{\emph{a}}[@webster1998] \\
            toluene&5.25&&\\
            \hline
        \end{tabular}
        \begin{tablenotes}[para,flushleft]
            \textsuperscript{\emph{a}}~Calculated critical diameter
        \end{tablenotes}
    \end{threeparttable}
    \label{stbl:test_molecules1}
\end{table}
\clearpage


# Application to research


Enzymes are a class of protein whose function is to catalyze the biochemical transformation of a substrate to a product.
Improved durability is a crucial step toward the commercial application of enzymes as industrial catalysts.[@schmid2001]
Through encapsulation in metal--organic frameworks (MOFs) matrices, enzymes are able to retain their activity in harsh conditions (e.g. elevated temperatures or proteolytic media) due to the protection afforded by the MOF matrix.[@liang2015]
Additionally, the MOF matrix is expected to afford size-selective transport of substrates to the active site of an enzyme via its pore network.
The cost of purified enzymes and the time required to prepare and screen new enzymes is significant.
Hence, accurate predictions of which reactions may be possible inside a given material before doing costly experiments are valuable.
Therefore, we implemented `mol-ellipsize` to enable efficient screening of large numbers of reactants and products of enzymatic reactions for their ability to diffuse within MOFs, based on their molecular size.
We used the KEGG database[@kanehisa2000; @kanehisa2016; @kanehisa2017] of $\sim$\num{12000} curated enzymatic reactions to compile a list of possible reactions with known molecular components.
From $\sim$\num{12000} KEGG reactions we extracted a chemical space of 5640 unique molecules (that were readable by RDKit with a mass limit of \SI{500}{\gram\per\mole}; a mass limit was applied because larger molecules are not expected to diffuse through typical MOFs).
Figure~\ref{fgr:chem_space} shows the calculated minimum intermediate diameter, $d$, of these molecules as a function of the number of heavy atoms in the molecule.
Unsurprisingly, there is a correlation between size and the number of heavy atoms, which begins to plateau as larger molecules tend to be more flexible and are, therefore, able to adopt conformations with smaller intermediate diameters than would be predicted from a naive consideration of their chemical formulae.
Figure~\ref{fgr:chem_space} clearly shows that while the majority of the chemical space is too large to diffuse through ZIF-8 (a prototypical porous material), there are still a significant number of molecules ($\sim$1000) with $d<\SI{6.6}{\angstrom}$ (the approximate diffusion threshold for ZIF-8)[@zhang2012; @zhang2013; @verploegh2015; @diestel2012; @zhu2016], which should be able to diffuse through ZIF-8.
Therefore, we would expect that novel candidate enzyme reactions can be found.
Furthermore, the generalisability of `mol-ellipsize} allows this approach to expand that chemical space through synthetic modifications to the framework or reaction components.

\begin{figure}[ht!]
\centering
\includegraphics{joss_chemspace}
\caption{
Distribution of the minimum intermediate diameter $d$ calculated by `mol-ellipsize} as a function of the number of heavy atoms in each molecule of all collected molecules.
The shaded region in indicates the approximate range for the threshold for diffusion through ZIF-8 from the literature.
}
\label{fgr:chem_space}
\end{figure}


All code used to extract molecules from KEGG is available at \url{https://github.com/andrewtarzia/enzyme_screen.git}.









# Acknowledgements

A.~T. was supported by an Australian Government RTP Scholarship and a PhD top-up scholarship from CSIRO Division of Materials Science and Engineering. A.~T. would like to thank Dr~Jesse Teo, Dr~Natasha Maddigan, Dr~Weibin Liang, Dr~Justin Spence, A/Prof~Stephen Bell, Dr~Austin Mroz and Dr~Kim E. Jelfs for discussions about this work.

# References






