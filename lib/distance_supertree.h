/* 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file distance_supertree.h
 *  \brief functions for carrying "matrix representation with distance" methods, including recent extensions for
 *  paralogous gene families
 *
 */

#ifndef _distance_supertree_h_
#define _distance_supertree_h_

#include <biomcmc.h> 

char_vector get_species_names_from_newick_space (newick_space g_nwk, char_vector spnames, bool check_spnames);
char_vector assume_species_names_from_newick_space (newick_space g_nwk);
newick_space find_matrix_distance_species_tree (newick_space g_nwk, char_vector spnames, double tolerance, bool check_spnames, bool remove_reorder_when_check_spnames, bool fast);

#endif
