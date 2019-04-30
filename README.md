# super\_distance (distance-based species supertrees)
[![Build Status](https://travis-ci.org/quadram-institute-bioscience/super_distance.svg?branch=master)](https://travis-ci.org/quadram-institute-bioscience/super_distance)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/quadram-institute-bioscience/super_distance/blob/master/LICENSE)
[![Docker Pulls](https://img.shields.io/docker/pulls/leomrtns/super_distance.svg)](https://hub.docker.com/r/leomrtns/super_distance)


This software implements the most common supertree methods, with emphasis on whole gene families (i.e. gene trees that
may contain paralogs) for species tree inference. The idea is to use a common framework for different algorithms, although not 
all methods are implemented (e.g.  [Bayesian supertrees](https://bitbucket.org/leomrtns/guenomu/) are absent).

The software is fast enough such that it can be added into workflows, and given a large set of input trees (that we call *gene
trees* or *gene families*) it produces a small set of output trees (*species trees*) that tries to summarise the
information in the input trees. 

### note for beta testers:
This software is still not complete, probably not ready for public consumption; the main limitations (specially when the program doesn't
behave like described in this readme) are marked with the symbol &#x26D4;.
Currently it assumes newick files, and uses only the distance-based estimation. 

## Installation
### Conda
[![Anaconda-Server Badge](https://anaconda.org/bioconda/super_distance/badges/platforms.svg)](https://anaconda.org/bioconda/super_distance)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/super_distance/badges/latest_release_date.svg)](https://anaconda.org/bioconda/super_distance)
After you install [miniconda](https://conda.io/en/latest/miniconda.html), simply run
```[bash]
conda install -c bioconda super_distance
```

### From source
The instalation uses the autotools build system for compilation, and relies on the
[biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) library, which can be downloaded
recursivelly:
```[bash]
/home/simpson/$ git clone --recursive git@github.com:quadram-institute-bioscience/super_distance.git
/home/simpson/$ mkdir build
/home/simpson/$ cd build
/home/simpson/$ ../super_distance/configure --prefix=/usr/local
/home/simpson/$ make; sudo make install
```
As seen above, it is usually good idea to compile the code on a dedicated clean directory (`build`, in the example). 
The example above will install the `libsuper_distance` globally, in `/usr/local/lib`. 
If you don't have administrative (*sudo*) priviledges you can chose a local directory, by replacing the two last lines
with:
```[bash]
/home/simpson/$ ../super_distance-master/configure --prefix=/home/simpson/local
/home/simpson/$ make; make install
```
You can then run a battery of tests with
```[bash]
/home/simpson/$ make check
```
If you download the zip instead of git-cloning you will miss the the biomcmc-lib library, which is a submodule. In this
case please [download it](https://github.com/quadram-institute-bioscience/biomcmc-lib) and unzip it below `super_distance-master/`.

### Docker
After installing [Docker](https://www.docker.com/), you can install a docker container for `super_distance` with:
```[bash]
docker pull leomrtns/super_distance
```
And to use it you can run something like  
```[bash]
docker run --rm -it -v /path/to/data:/data leomrtns/super_distance sh -c 'super_distance -s /data/species_names.txt /data/gene*.tre'
```
Notice that the command we invoke is actually `sh`, to be able to use shell expansion; you are free to call
`super_distance` directly but in this case you must write all file names (thanks to [Andrea](https://github.com/telatin) for the trick!).
You may also prefer to first go to the working directory (`/path/to/data` in our example) and then run everything from
there &mdash; remember that docker can only access files mounted with `-v`:
```[bash]
docker run --rm -it -v /path/to/data:/data leomrtns/super_distance sh -c 'cd /data && super_distance -s species_names.txt gene*.tre'
```

## Usage 
The program needs a file with the species names (see below) and a list of the gene trees, in newick format. You can see
the options by running
```[bash]
/home/simpson/$ super_distance -h      ## OR
/home/simpson/$ super_distance --help  ## same as above
```
As seen above, all options have a short (one character) and a long version. Currently the available options are:
- **--epsilon (-e)** This is the minimum branch length on internodal distances; values smaller than this are considered to
  be multifurcations
- **--species (-s)** file name with the list of species names
- **--output (-o)** output filename with resulting supertrees
- **--mode (-m)** which algorithms should be used, with a one-letter code for each. This is ever-changing at the moment,
  only the distance-based algortihms are certain. 
- all remaining arguments are assumed to be the names of gene files. 
The list above may be incomplete as the software is being developed. Please run `super_distance -n` for an up-to-date
description. 

### Mapping from genes to species
The input trees don't need to have information on all species, which is the classic supertree setting. 
They can also have the same species represented more than once, as when we have whole gene families with paralogs, 
and/or several samples from same species as in population genomics data sets. 
These trees with the same label for several leaves are called "multi-labelled trees", or simply *mul-trees*.
We use this term when we want to emphasise the distinction from classic supertrees approaches (where the objective was
to create a tree on the full set of taxa, from trees on subsets of it).
super\_distance works as expected on the classic setting, by the way. 

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
In the newick file it is valid to have several leaves with the same name, e.g. the species name, although most other software 
won't allow it.
Nexus files need unique taxon names, since the nexus format may need to map a sequence to a tree leaf.
(&#x26D4; currently only newick files are accepted, altough we alread have the equivalent functions for nexus trees).
`Super_distance` does not respect, however, spaces within a gene leaf or species names, despite these
agreeing with the [formal newick specification](http://evolution.genetics.washington.edu/phylip/newick_doc.html).
Remember that the list of species names will define the leaves of the output trees, and thus the program may work even
in the presence of spaces (since it removes them from all input files).
In the future, and if it bothers enough people, we may implement automatic inference of species names. 
It is worh mentioning again that the software works equally well in the absence of mul-trees, but the file with species
names must still be provided. 

## Algorithms
Several distance-based and one bipartition-based supertree methods are being implemented, but only the distance based
methods are functional. 
Multifurcating trees are allowed, since the polytomies are transformed into dicotomies of length zero. 

### Distance-based, or MRD supertrees
These methods are sometimes called "matrix representation with distances" (MRD), specially in the classic supertree
context (where we don't have mul-trees), and are a generalisation of the ASTRID, NJst, STAR, and a few others. 
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
or the pairwise species matrix. 
To avoid a [Buridan's donkey situation](https://en.wikipedia.org/wiki/Buridan%27s_ass), we've implemented all possible
combinations, with a few caveats: 
  
  1. We only use the average, and not the minimum, between loci. Within a locus (gene) we can use both. 
  2. We implemented both UPGMA and single-linkage clustering, besides the bioNJ implemenentation of the
     Neighbour-Joining algorithm. 
  3. When we rescale the gene trees, we scale back the final pairwise distance matrix, before the clustering step. This
     final scaling is based on the average over all genes, such that all supertrees should have easily interpretable lengths
     (except maybe for the internodal distances). 
  4. It is not uncommon to have a lot of missing information, for instance when two species are never seen together in
     the same gene. In this case we estimate their pairwise distance from species in common using the [ultrametric
     approach](http://dx.doi.org/10.1093/bioinformatics/bth211).
As usual, some methods/combinations will make more sense than others. 

&#x26D4; Currently we report a list of supertrees without any explanation about the method, but this may change (if
we allow the user to set the models, although I prefer to infer all since they're fast)

### Bipartition-based, or MRP supertrees
These are the classic supertree approaches, also known as "matrix representation with parsimony" (MRP) since the
maximum parsimony tree is inferred from the bipartition patterns.
However we extended it to work with mul-trees, by looking at the species represented at both ends of each bipartition
(&#x26D4; *unfinished*, right now it works correctly only with ortholog sets, and is not offered to user).

The set of gene trees will generate a binary matrix where each row (sample) is a species and each column (dimension) is a 
bipartition.
For any given species tree the parsimony score can be calculated from this matrix. where missing data is coded
accordingly, and a initial state can be found through a clustering algorithm.
This is a preferred algorithm for sparse data sets (i.e. gene trees don't have information about all species). 
Several modifications of the basic algorithm are possible, as the maximum compatibility supertree or the quartet
supertree (although we don't have plans of implementing the latter). 

### Reconciliation, or median supertrees
These are supertrees that try to minimise the distance from the set of input trees.
They can be with respect to a particular tree-to-tree distance, or to a set of distances &mdash; in which case no single
supertree will be a global optimum.
Right now they neglect branch lengths, but work seamlessly with mul-trees. 
The quartet supertree would be implemented here (as the ASTRAL algorithm). 
(&#x26D4; *unfinished*, currently the user cannot chose this model).

### Consensus trees
If all the input trees share the same leaf set (i.e. the same species, with no missing data), then we can estimate the
consensus trees.
In all supertree methods above we assume one sample tree per gene family, but here we can create weightings per gene
tree file &mdash; nexus files are particularly suited for large collections since they have a translation table for leaf
names, and allow for compact distributions of topologies (as in the `.trprobs` files of MrBayes, for instance)
(&#x26D4; *unfinished*, currently the user cannot chose this mode; I would offer it as another program).

## License 
Copyright (C) 2019-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

super\_distance is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

This software borrows many functions from the [guenomu software](https://bitbucket.org/leomrtns/guenomu/) for phylogenomic species tree inference and 
also extends functionality from [genefam-dist library](https://github.com/leomrtns/genefam-dist). It relies on the
low-level library [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib), which is defined as a submodule &mdash; so don't forget to git recursively. 
