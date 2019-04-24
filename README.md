# super_sptree (species supertrees)

This software implements the most common supertree methods, with emphasis on whole gene families (i.e. gene trees that
may contain paralogs). The idea is to use a common framework for different algorithms, although not all methods are
implemented (e.g.  [Bayesian supertrees](https://bitbucket.org/leomrtns/guenomu/) are absent).

The software is fast enough such that it can be added into workflows, and given a large set of input trees (that we call *gene
trees* or *gene families*) it produces a small set of output trees (*species trees*) that tries to summarise the
information in the input trees. 

### Mapping from genes to species
The input trees don't need to have information on all species, which is the classic supertree setting. 
They can also have the same species represented more than once, as when we have whole gene families with paralogs, 
and/or several samples from same species as in population genomics data sets. 
These trees with the same label for several leaves are called "multi-labelled trees", or simply *mul-trees*.
We use this term when we want to emphasise the distinction from classic supertrees approaches (where the objective was
to create a tree on the full set of taxa, from trees on subsets of it).

Therefore, besides the input gene trees the program will request a file with a list of species names, which will provide
a mapping between leaves from the gene trees and leaves from the species tree. 
The program tries to find the species name associated to each gene leaf by string matching, which means the gene tree
leaves must contain the species names. For instance if the list of species names is

```[bash]
Neisseria_elongata
Neisseria_gonorrhoeae
Neisseria_lactamica
Neisseria_meningitidis
Pseudomonas_aeruginosa
Pseudomonas_brassicacearum
Pseudomonas_chlororaphis
```

Then the gene leaves would be mapped as follow:
```[bash]
|      gene leaf names           | species it will be mapped to  |
|--------------------------------|-------------------------------|
| Neisseria_elongata_001         | Neisseria_elongata            |
| Neisseria_elongata_002         | Neisseria_elongata            |
| COG001_Neisseria_gonorrhoeae   | Neisseria_gonorrhoeae         |  
| Neisseria_gonorrhoeae_COI2     | Neisseria_gonorrhoeae         |   
| Pseudomonas_chlororaphis       | Pseudomonas_chlororaphis      |   
| Pseudomonas_brassicacearum     | Pseudomonas_brassicacearum    |   
| Pseudomonas_brassica           | <<UNKNOWN>>                   | 
```
It is valid to have several leaves with the same name, e.g. the species name, although most other software won't allow it.
`Super_sptree` does not respect, however, spaces within a gene leaf or species names, despite
agreeing with the [formal newick specification](http://evolution.genetics.washington.edu/phylip/newick_doc.html).
Remember that the list of species names will define the leaves of the output trees.
In the future, and if it bothers enough people, we may implement automatic inference of species names. 
It is worh mentioning that the software works equally well in the absence of mul-trees, and some extra functionality is
available in such cases.

## Algorithms

Currently several distance-based and one bipartition-based supertree methods are implemented.

### Distance-based
These methods are sometimes called "matrix representation with distances" (MRD), specially in the classic supertree
context (where we don't have mul-trees). and are a generalisation of the ASTRID,
NJst, STAR, and a few others. 
What they all have in common is that 

  1. For each gene tree they create a matrix with 'patristic' distance between leaves
  2. If there are several leaves from same species (i.e. the gene tree is a mul-tree) then the average or the minimum over all
     possible patristic distances is taken as the species-wise distance between pairs
  3. Once they have one pairwise distance matrix per gene, they merge them into one matrix by taking the average or the
     minimum across genes. 
  4. This overall distance matrix is then used by a clustering algorithm to estimate the species tree.

Some methods use rooted while others use unrooted trees; some use branch lengths while
others use just the internodal distances; some resolve conflicts by taking the average or the minimum distances, 
within or between gene trees; some use UPGMA and some use NJ to estimate the species tree; some rescale branch lengths
or the pairwise species tree. 
To avoid a [Buridan's donkey situation](https://en.wikipedia.org/wiki/Buridan%27s_ass), we've implemented all possible
combinations, with a few caveats: 
  
  1. We only use the average, and not the minimum, between loci. Within a locus (gene) we can use both. 
  2. We implemented both UPGMA and single-linkage clustering, besides the bioNJ implemenentation of the
     Neighbour-Joining algorithm. 
  3. We always scale the final pairwise distance matrix, before the clustering step. Branch lengths in original scale
     can be recovered from subsets of gene trees (*not implemented yet*). 
  4. It is not uncommon to have a lot of missing information, for instance when two species are never seen together in
     the same gene. In this case we estimate their pairwise distance from species in common (*not implemented yet*).

As usual, some methods/combinations will make more sense than others.

## Installation
### From source
The instalation uses the autotools build system for compilation, and relies on the
[biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) library, which can be downloaded
recursivelly:
```[bash]
/home/simpson/$ git clone --recursive git@github.com:quadram-institute-bioscience/super_sptree.git
/home/simpson/$ mkdir build
/home/simpson/$ cd build
/home/simpson/$ ../super_sptree-master/configure --prefix=/usr/local
/home/simpson/$ make; sudo make install
```

As seen above, it is usually good idea to compile the code on a dedicated clean directory (`build`, in the example). 
The example above will install the `libsuper_sptree` globally, in `/usr/local`. 
If you don't have administrative (*sudo*) priviledges you can chose a local directory, by replacing the two last lines
with:
```[bash]
/home/simpson/$ ../super_sptree-master/configure --prefix=/home/simpson/local
/home/simpson/$ make; make install
```
You can then run a battery of tests with
```[bash]
/home/simpson/$ make check
```

If you download the zip instead of git-cloning you will miss the the biomcmc-lib library, which is a submodule. In this
case please [download it](https://github.com/quadram-institute-bioscience/biomcmc-lib) and unzip it below `super_sptree-master/`.

## License 
Copyright (C) 2019-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

super_sptree is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

This software borrows many functions from the [guenomu software](https://bitbucket.org/leomrtns/guenomu/) for phylogenomic species tree inference and 
also extends functionality from [genefam-dist library](https://github.com/leomrtns/genefam-dist). It relies on the
low-level library [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib), which is defined as a submodule &mdash; so don't forget to git recursively. 
