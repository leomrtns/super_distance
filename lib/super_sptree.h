/* 
 * This file is part of genefam-dist, a library for calculating distances between gene families and species trees. 
 * Copyright (C) 2016  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * genefam-dist is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file super_sptree.h 
 *  \brief biomcmc library interface to external programs, specific to super_sptree repo.
 *
 *  The idea is for super_sptree is to be general for several sofware, including treesignal but also parallel projects by
 *  leomrtns. This library started as a branching from the biomcmc library (from the guenomu software)
 */

#ifndef _super_sptree_h_
#define _super_sptree_h_

#include <biomcmc.h> 
#include "topology_space.h"
#include "topology_mrca.h"
#include "topology_build.h"

#ifdef THESE_ARE_COMMENTS
#include "topology_splitset.h" // called by topology_space.h
#endif // of THESE_ARE_COMMENTS

/*! \brief given a gene tree and a group of species trees, both in newick format, return the spectrum of unnormalized distances */
int genefam_module_treesignal_fromtrees (const char *gtree_str, const char *splist_str, double **output_distances);
/*! \brief given a gene tree and a group of species trees, both in newick format, return the spectrum of rescaled distances */
int genefam_module_treesignal_fromtrees_rescale (const char *gtree_str, const char *splist_str, double **output_distances);
/*! \brief given a gene tree and a group of species trees, both in newick format, return frequencies of random sptrees with distances smaller than group */
int genefam_module_treesignal_fromtrees_pvalue (const char *gtree_str, const char *splist_str, int n_reps, double **output_distances);
/*! \brief given set of trees, return trees with SPR neighbours (as newick string). Useful for generating reference sptrees */
char* genefam_module_randomise_trees_with_spr (const char *splist_str, int n_copies, int n_spr);
/*! \brief generate a set of trees (as newick string), all separated from previous one by a number of SPR moves */
char* genefam_module_generate_spr_trees (int n_leaves, int n_iter, int n_spr);

#endif
