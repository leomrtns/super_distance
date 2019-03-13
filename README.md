# super_sptree (super species tree)

This software is written in C, with functions to calculate distances between trees and to find no-so-bad supertrees (species trees, usually). 
It borrows many functions from the [guenomu software](https://bitbucket.org/leomrtns/guenomu/) for phylogenomic species tree inference and 
also from the [genefam-dist library](https://github.com/leomrtns/genefam-dist).

The code relies on [biomcmc-lib](https://github.com/leomrtns/biomcmc-lib) which is defined as a submodule --- so don't forget to git recursively. 

## Installation

The instalation uses the autotools build system.  

```
/home/simpson/$ git clone --recursive git@github.com:leomrtns/super_sptree.git
/home/simpson/$ mkdir build
/home/simpson/$ cd build
/home/simpson/$ ../super_sptree-master/configure --prefix=/usr/local
/home/simpson/$ make; sudo make install
```

As seen above, it is usually good idea to compile the code on a dedicated clean directory (`build`, in the example). For
some reason many other software using autotools don't work properly, forcing you to run `configure` inside the
repository directory...

The example above will install the `libsuper_sptree` globally, in `/usr/local`. If you don't have administrative priviledges you can
chose a local directory (and drop the "sudo" command).

BTW if you download the zip instead of cloning you may miss the the biomcmc-lib library, which is a submodule.

## License 
Copyright (C) 2019-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

super_sptree is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

