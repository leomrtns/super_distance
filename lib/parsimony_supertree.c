#include "parsimony_supertree.h" 

void find_upgma_mrp_tree_and_add_to_newick_space (double *dists, char_vector spnames, newick_space nwk, int tree_method);

newick_space
find_upgma_mrp_species_tree (newick_space g_nwk, char_vector spnames, bool check_spnames, bool remove_reorder_when_check_spnames)
{
  int i, *sp_idx_in_gene = NULL;
  double *d_w = NULL, *d_u = NULL;
  char_vector species_names;
  binary_parsimony pars;
  newick_space species_nwk = new_newick_space ();

  /* 1. remove species absent from all genes */
  if (check_spnames) species_names = get_species_names_from_newick_space (g_nwk, spnames, remove_reorder_when_check_spnames);
  else {
    species_names = spnames;
    spnames->ref_counter++;
  }
  pars = new_binary_parsimony (species_names->nstrings); /* no sptree information, just a matrix from genetrees */

  /* 2. update mrp columns with bipartitions from trees */
  for (i=0; i < g_nwk->ntrees; i++) {
    sp_idx_in_gene = (int *) biomcmc_realloc ((int*) sp_idx_in_gene, g_nwk->t[i]->nleaves * sizeof (int));
    index_species_gene_char_vectors (species_names, g_nwk->t[i]->taxlabel, sp_idx_in_gene, NULL);
    update_binary_parsimony_from_topology (pars, g_nwk->t[i], sp_idx_in_gene);
  }
  /* 3. create square distance_matrix */
  d_w = (double *) biomcmc_malloc ((species_names->nstrings * (species_names->nstrings - 1)/2) * sizeof (double));
  d_u = (double *) biomcmc_malloc ((species_names->nstrings * (species_names->nstrings - 1)/2) * sizeof (double));

  pairwise_distances_from_binary_parsimony_datamatrix (pars->external, d_w, d_u);

  /* 3. find upgma and bionj trees, for both unweighted and weighted distance matrices */
  for (i = 0; i < 3; i++) {
    find_upgma_mrp_tree_and_add_to_newick_space (d_w, species_names, species_nwk, i);
    find_upgma_mrp_tree_and_add_to_newick_space (d_u, species_names, species_nwk, i);
  }

  if (d_w) free (d_w);
  if (d_u) free (d_u);
  if (sp_idx_in_gene) free (sp_idx_in_gene);
  del_char_vector (species_names);
  return species_nwk;
}

void
find_upgma_mrp_tree_and_add_to_newick_space (double *dists, char_vector spnames, newick_space nwk, int tree_method)
{
  int i,j;
  topology maxtree;
  distance_matrix square;

  maxtree = new_topology (spnames->nstrings);
  maxtree->taxlabel = spnames; spnames->ref_counter++; /* sptaxa is pointed here and at the calling function */
  square = new_distance_matrix (spnames->nstrings);
  for (j=1; j < spnames->nstrings; j++) for (i = 0; i < j; i++) square->d[i][j] = dists[ ((j * (j-1)) / 2 + i) ];

  if (tree_method == 0) bionj_from_distance_matrix (maxtree, square); 
  else if (tree_method == 1) upgma_from_distance_matrix (maxtree, square, false); // false -> upgma, true -> single linkage 
  else upgma_from_distance_matrix (maxtree, square, true); // false -> upgma, true -> single linkage 
  update_newick_space_from_topology (nwk, maxtree);
  del_distance_matrix (square);
  return;
}

