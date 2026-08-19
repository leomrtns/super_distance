#include "distance_supertree.h" 

void find_maxtree_and_add_to_newick_space (spdist_matrix dist, char_vector spnames, newick_space nwk, int tree_method, bool use_within_gf_means);
void add_tree_fitted_distances (spdist_matrix dist, topology tree, newick_space nwk);
static bool newick_space_has_paralogs (newick_space g_nwk, char_vector spnames);

char_vector
get_species_names_from_newick_space (newick_space g_nwk, char_vector spnames, bool remove_reorder)
{
  int i, *sp_count, *valid, n_valid = 0;
  char_vector vec;

  if (remove_reorder) { // will change original spnames #char_vector 
    char_vector_remove_duplicate_strings (spnames); /* duplicate names */
    char_vector_reorder_by_size_or_lexicographically (spnames, false, NULL); // false/true -> by size/lexico
  }

  sp_count = (int *) biomcmc_malloc (spnames->nstrings * sizeof (int));
  valid = (int *) biomcmc_malloc (spnames->nstrings * sizeof (int));
  for (i=0; i < spnames->nstrings; i++) valid[i] = sp_count[i] = 0;
  for (i=0; i < g_nwk->ntrees; i++) update_species_count_from_gene_char_vector (spnames, g_nwk->t[i]->taxlabel, sp_count);
  for (i=0; i < spnames->nstrings; i++) if (sp_count[i] > 0) valid[n_valid++] = i;
  if (n_valid == spnames->nstrings) {
    free (sp_count);
    free (valid);
    spnames->ref_counter++;
    return spnames;
  }
  vec = new_char_vector_from_valid_strings_char_vector (spnames, valid, n_valid);
  char_vector_reorder_by_size_or_lexicographically (vec, false, NULL); // false/true -> by size/lexico
  if (sp_count) free (sp_count);
  if (valid) free (valid);
  return vec;
}

char_vector
assume_species_names_from_newick_space (newick_space g_nwk)
{
  int i, j, nleaves = 0;
  char_vector vec;
  hashtable species_hash;
  /* collect the union of leaf names across all trees (assumes orthology) */
  for (i = 0; i < g_nwk->ntrees; i++) nleaves += g_nwk->t[i]->taxlabel->nstrings;
  vec = new_char_vector (nleaves);
  species_hash = new_hashtable (nleaves);
  for (i = 0; i < g_nwk->ntrees; i++)
    for (j = 0; j < g_nwk->t[i]->taxlabel->nstrings; j++) if (lookup_hashtable (species_hash, g_nwk->t[i]->taxlabel->string[j]) < 0) {
      char_vector_add_string (vec, g_nwk->t[i]->taxlabel->string[j]);
      insert_hashtable (species_hash, g_nwk->t[i]->taxlabel->string[j], vec->next_avail - 1);
    }
  del_hashtable (species_hash);
  char_vector_remove_empty_strings (vec);
  char_vector_reorder_by_size_or_lexicographically (vec, false, NULL); // false/true -> by size/lexico
  return vec;
}

static bool
newick_space_has_paralogs (newick_space g_nwk, char_vector spnames)
{
  int i, j, *sp_idx_in_gene, *last_tree_seen;
  bool has_paralogs = false;

  last_tree_seen = (int*) biomcmc_malloc (spnames->nstrings * sizeof (int));
  for (i = 0; i < spnames->nstrings; i++) last_tree_seen[i] = -1;
  for (i = 0; (i < g_nwk->ntrees) && (!has_paralogs); i++) {
    sp_idx_in_gene = (int*) biomcmc_malloc (g_nwk->t[i]->nleaves * sizeof (int));
    index_species_gene_char_vectors (spnames, g_nwk->t[i]->taxlabel, sp_idx_in_gene, NULL);
    for (j = 0; (j < g_nwk->t[i]->nleaves) && (!has_paralogs); j++) {
      if (last_tree_seen[sp_idx_in_gene[j]] == i) has_paralogs = true;
      else last_tree_seen[sp_idx_in_gene[j]] = i;
    }
    free (sp_idx_in_gene);
  }
  free (last_tree_seen);
  return has_paralogs;
}

newick_space
find_matrix_distance_species_tree (newick_space g_nwk, char_vector spnames, double tolerance, bool check_spnames, bool remove_reorder_when_check_spnames, bool fast, bool assume_orthology)
{
  int i, j, k, n_pairs, n_dists = 6, n_trees_to_rescale, *sp_idx_in_gene = NULL;
  double **dist, *this_scaling, *scaling;
  bool has_paralogs;
  char_vector species_names;
  spdist_matrix *dm_glob, dm_local; 
  newick_space species_nwk = new_newick_space ();

  if (fast) n_dists = 2; // only nodal and average
  /* 1. remove species absent from all genes */
  if (check_spnames) species_names = get_species_names_from_newick_space (g_nwk, spnames, remove_reorder_when_check_spnames);
  else {
    species_names = spnames;
    spnames->ref_counter++;
  }
  has_paralogs = assume_orthology ? false : newick_space_has_paralogs (g_nwk, species_names);
  /* 1.5 create structures, remembering that dm_glob have only _means_ across loci */
  dm_glob = (spdist_matrix*) biomcmc_malloc (n_dists * sizeof (spdist_matrix));
  for (j=0; j < n_dists; j++) dm_glob[j] = new_spdist_matrix (species_names->nstrings);
  for (j=0; j < n_dists; j++) zero_all_spdist_matrix (dm_glob[j]); // zero min[] since we'll calc the average in the end  
  dist = (double**) biomcmc_malloc (n_dists * sizeof (double*)); 
  for (j=0; j < n_dists; j++) dist[j] = NULL;
  this_scaling = (double*) biomcmc_malloc (n_dists * sizeof (double)); 
  scaling = (double*) biomcmc_malloc (n_dists * sizeof (double)); 
  for (j=0; j < n_dists; j++) scaling[j] = this_scaling[j] = 0.; 
  dm_local  = new_spdist_matrix (species_names->nstrings);

  /* 2. update species distance matrices */
  for (i=0; i < g_nwk->ntrees; i++) {
    n_pairs = (g_nwk->t[i]->nleaves * (g_nwk->t[i]->nleaves -1))/2;
    sp_idx_in_gene = (int*) biomcmc_malloc (g_nwk->t[i]->nleaves * sizeof (int));
    index_species_gene_char_vectors (species_names, g_nwk->t[i]->taxlabel, sp_idx_in_gene, NULL);

    for (j=0; j < n_dists; j++) dist[j] = (double*) biomcmc_realloc ((double*) dist[j], n_pairs * sizeof (double));
    patristic_distances_from_topology_to_vectors (g_nwk->t[i], dist, this_scaling, n_dists, tolerance); 
    for (j=0; j < n_dists; j++) {
      fill_spdistmatrix_from_gene_dist_vector (dm_local, dist[j], g_nwk->t[i]->nleaves, sp_idx_in_gene);
      update_spdistmatrix_from_spdistmatrix (dm_glob[j], dm_local);
      scaling[j] += this_scaling[j];
    }
    if (sp_idx_in_gene) free (sp_idx_in_gene);
  }
  for (j=0; j < n_dists; j++) {
    scaling[j] /= (double)(g_nwk->ntrees);
    finalise_spdist_matrix_with_rescaling (dm_glob[j], scaling[j]);
  }

  if (dm_glob[0]->n_missing) fprintf (stderr, "OBS: %d species pair combinations never appear on same gene family\n", dm_glob[0]->n_missing);
  /* 3. find upgma and bionj trees, for both unweighted and weighted distance matrices */

  if (fast) {
    for (j=0; j < n_dists; j++)
      find_maxtree_and_add_to_newick_space (dm_glob[j], species_names, species_nwk, 1, true); // 1->upgma; true->means within locus
    estimate_topology_branch_lengths_from_distances (species_nwk->t[0], dm_glob[1]->mean);// average branch lengths applied to nodal distance tree
  }
  else for (j=0; j < n_dists; j++) for (i = 0; i < 3; i++) { 
    find_maxtree_and_add_to_newick_space (dm_glob[j], species_names, species_nwk, i, true);
    if (has_paralogs) find_maxtree_and_add_to_newick_space (dm_glob[j], species_names, species_nwk, i, false);
    if (j == 0) { // trees from nodal distances are rescaled from average branch lengths
      n_trees_to_rescale = has_paralogs ? 2 : 1;
      for (k = n_trees_to_rescale; k > 0; k--)
        estimate_topology_branch_lengths_from_distances (species_nwk->t[species_nwk->ntrees-k], dm_glob[1]->mean);
    }
  }

  if (dist) {
    for (j = n_dists-1; j >=0; j--) if (dist[j]) free (dist[j]);
    free (dist);
  }
  for (j = n_dists-1; j >=0; j--) del_spdist_matrix (dm_glob[j]);
  if (this_scaling) free (this_scaling);
  if (scaling) free (scaling);
  if (dm_glob) free (dm_glob);
  del_spdist_matrix (dm_local);
  del_char_vector (species_names);
  return species_nwk;
}

void
find_maxtree_and_add_to_newick_space (spdist_matrix dist, char_vector spnames, newick_space nwk, int tree_method, bool use_within_gf_means)
{
  topology maxtree;
  distance_matrix square = new_distance_matrix (spnames->nstrings);
  copy_spdist_matrix_to_distance_matrix_upper (dist, square, use_within_gf_means);
  maxtree = new_topology (spnames->nstrings);
  maxtree->taxlabel = spnames; spnames->ref_counter++; /* sptaxa is pointed here and at the calling function */

  if (tree_method == 0) bionj_from_distance_matrix (maxtree, square); 
  else if (tree_method == 1) upgma_from_distance_matrix (maxtree, square, false); // false -> upgma, true -> single linkage 
  else upgma_from_distance_matrix (maxtree, square, true); // false -> upgma, true -> single linkage 
  update_newick_space_from_topology (nwk, maxtree);
  del_distance_matrix (square);
  return;
}

void
add_tree_fitted_distances (spdist_matrix dist, topology tree, newick_space nwk)
{
  topology this = new_topology (tree->nleaves);
  copy_topology_from_topology (this, tree);
  estimate_topology_branch_lengths_from_distances (this, dist->mean);
  update_newick_space_from_topology (nwk, this);
}
