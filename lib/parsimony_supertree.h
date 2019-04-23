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

/*! \file parsimony_supertree.h
 *  \brief functions for carrying "matrix representation with parsimony" methods, including experimental extensions to
 *  paralogous gene families
 *
 */
#ifndef _parsimony_supertree_h_
#define _parsimony_supertree_h_

#include "distance_supertree.h" 

/*! \brief initial MRP tree is found by nj/upgma on binary (bipartition matrices). #check_spnames true if
 * get_species_names_from_newick_space() called before but now we are on subset. In this case
 * #remove_reorder_when_check_spnames can be false */
newick_space find_upgma_mrp_species_tree (newick_space g_nwk, char_vector spnames, bool check_spnames, bool remove_reorder_when_check_spnames);

#endif
