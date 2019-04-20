#include <super_sptree.h> 
/* based on guenomu/maxtree but also implementing elements from guenomu in ML mode */

char_vector
get_species_names_from_newick_space (newick_space g_nwk, char_vector spnames)
{
  int i, *sp_count, *valid, n_valid = 0;
  char_vector vec;
  sp_count = (int *) biomcmc_malloc (spnames->nstrings * sizeof (int));
  for (i=0; i < spnames->nstrings; i++) sp_count[i] = 0;
  for (i=0; i < g_nwk->ntrees; i++) update_species_count_from_gene_char_vector (spnames, g_nwk->t[i]->taxlabel, sp_count);
  for (i=0; i < spnames->nstrings; i++) if (sp_count[i] > 0) valid[n_valid++] = i;
  if (n_valid == spnames->nstrings) return spnames;
  vec = new_char_vector_from_valid_strings_char_vector (spnames, int *valid, int n_valid);
  char_vector_reorder_by_size_or_lexicographically (vec, false, NULL); // false/true -> by size/lexico
  return vec;
}

newick_space
find_matrix_distance_species_tree (newick_space g_nwk, char_vector spnames)
{
  char_vector species_names;
  species_names = get_species_names_from_newick_space (g_nwk, spnames);
  if (species_names == spnames) spnames->ref_counter++;



}
 
int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i,j,k, n_samples = 200, n_cycles = 2, n_iter_per_cycle = 16, n_output_trees = 2, n_ratchet = 128, n_proposal = 4, fraction_gfs = 512;
  char_vector species, genefiles;
  topology_space tsp;
  gene_sptrees gs;
  FILE *stream;
  char *s;

  time0 = clock ();
  print_welcome (argv[0]);
  if (argc < 3) print_usage (argv[0]);
  biomcmc_random_number_init(0);

  species = new_char_vector_from_file (argv[1]);
  genefiles = new_char_vector (1); // updated by maxtrees_from_subsamples()
  char_vector_remove_duplicate_strings (species); /* duplicate names */
  char_vector_longer_first_order (species); // order species names from longest to shortest (speed up gene/spnames comparison) 
  tsp = maxtrees_from_subsamples (species, argv + 2, argc - 2, genefiles); 
  save_topology_space_to_trprobs_file (tsp, "patristic.tre", 1.);

	return EXIT_SUCCESS;
}

topology_space
maxtrees_from_subsamples (char_vector sptaxa, char **arg_filename, int arg_nfiles, char_vector gfilename)
{
  int i, *idx_gene_to_sp = NULL, valid_species_size = 0;
  topology_space tsp, genetre;
  pooled_matrix pool;
  /* order species names from longer names to shorter (so that EColi is searched only after EColiII, e.g.) */
  empfreq ef = new_empfreq_sort_decreasing (sptaxa->nchars, sptaxa->nstrings, 1); /* 1=size_t (0=char, 2=int)*/

  pool = new_pooled_matrix ((int)(arg_nfiles/10), sptaxa->nstrings);

  for (i = 0; i < arg_nfiles; i++) { /* scan through gene files */
    genetre  = read_topology_space_from_file (arg_filename[i], NULL, false);// read i-th gene family; false -> neglects rooting
    idx_gene_to_sp = (int *) biomcmc_realloc ((int*) idx_gene_to_sp, genetre->distinct[0]->nleaves * sizeof (int));/* for each gene leaf, index of species */
    index_sptaxa_to_genetaxa (sptaxa, genetre->taxlabel, idx_gene_to_sp, ef);/* map species names to gene names and store into idx[] */
    valid_species_size = prepare_spdistmatrix_from_gene_species_map (pool->this_gene_spdist, idx_gene_to_sp, genetre->distinct[0]->nleaves);
    if (valid_species_size > 4) {
      char_vector_add_string (gfilename, arg_filename[i]);
      update_pooled_matrix_from_gene_tree (pool, genetre->distinct[0], idx_gene_to_sp);
    }
    del_topology_space (genetre); 
  }

  fprintf (stderr, " - Using %d valid gene family tree files from the input set of %d files ",gfilename->nstrings, arg_nfiles);
  fprintf (stderr, "(%d were too small or had less than 4 distinct species)\n", arg_nfiles - gfilename->nstrings);
  fprintf (stderr, " - Creating %d distance matrices with each gene family participating in %d random ones\n\n", pool->n_sptrees, pool->n_sets_per_gene);
  finalise_pooled_matrix (pool, sptaxa);

  tsp = new_topology_space ();
  tsp->taxlabel = sptaxa; sptaxa->ref_counter++; /* sptaxa is shared here as well */
  for (i =0; i < 3; i++) {
    find_maxtree_and_add_to_topol_space (pool->d_total, pool, tsp, i, true);  // mean length within gene family (locus)
    find_maxtree_and_add_to_topol_space (pool->d_total, pool, tsp, i, false); // min lengths within gene family
  }
  for (i = 0; i < pool->n_sptrees; i++) {
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, 0, true); // bionj/upgma/singlelinkage, means/mins
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, 1, false); // upgma on mins 
    find_maxtree_and_add_to_topol_space (pool->d[i], pool, tsp, 2, false); // single linkage on mins 
  }
  del_empfreq (ef);
  del_pooled_matrix (pool);
  if (idx_gene_to_sp) free (idx_gene_to_sp);
  return tsp;
}

/* pooled_matrix_struct functions */

pooled_matrix
new_pooled_matrix (int n_sets, int n_species)
{
  int i;
  pooled_matrix pool;

  pool = (pooled_matrix) biomcmc_malloc (sizeof (struct pooled_matrix_struct));

  pool->n_sptrees = n_sets;
  pool->n_species = n_species;
  pool->n_sets_per_gene = n_sets/4;
  pool->next = 0; // vector is shuffled whenever next arrives at last element

  if (pool->n_sptrees < 1) pool->n_sptrees = 1; 
  if (pool->n_sets_per_gene < 1) pool->n_sets_per_gene = 1;

  pool->this_gene_spdist = new_spdist_matrix (n_species);
  pool->d_total = new_spdist_matrix (n_species);
  pool->square_matrix = new_distance_matrix (n_species);
  pool->d = (spdist_matrix*) biomcmc_malloc (pool->n_sptrees * sizeof (spdist_matrix));
  for (i = 0; i < pool->n_sptrees; i++) {
    pool->d[i] = new_spdist_matrix(n_species);
    zero_all_spdist_matrix (pool->d[i]);
  }
  return pool;
}

void
del_pooled_matrix (pooled_matrix pool)
{
  int i;
  if (! pool) return;
  if (pool->d) {
    for (i = pool->n_sptrees-1; i >=0; i--) del_spdist_matrix (pool->d[i]);
    free (pool->d);
  }
  del_spdist_matrix (pool->this_gene_spdist);
  del_spdist_matrix (pool->d_total);
  del_distance_matrix (pool->square_matrix);
  free (pool);
  return;
}

int
next_element_from_pooled_matrix (pooled_matrix pool)
{
  if (++pool->next < pool->n_sptrees) return pool->next;

  int i, j;
  spdist_matrix pivot;
  for (i = pool->n_sptrees - 1; i > 0; i--) {  // Knuth shuffle
    j = biomcmc_rng_unif_int (i+1); // including i
    pivot = pool->d[j];
    pool->d[j] = pool->d[i];
    pool->d[i] = pivot;
  }
  pool->next = 0;
  return pool->next;
}

void
update_pooled_matrix_from_gene_tree (pooled_matrix pool, topology gene_topol, int *idx_gene_to_species)
{
  int i, j;
  distance_matrix genedist = new_distance_matrix_for_topology (gene_topol->nleaves);
  // here we can play with branch lengths (now it's NULL)
  fill_distance_matrix_from_topology (genedist, gene_topol, NULL, true); /* true=store in upper diagonal */
  fill_spdistmatrix_from_gene_dists (pool->this_gene_spdist, genedist, idx_gene_to_species, true);
  for (i = 0; i < pool->n_sets_per_gene; i++) {
    j = next_element_from_pooled_matrix (pool); // shuffles pool->d[] when necessary 
    update_spdistmatrix_from_spdistmatrix (pool->d[j], pool->this_gene_spdist);
  }
  update_spdistmatrix_from_spdistmatrix (pool->d_total, pool->this_gene_spdist);
  del_distance_matrix (genedist);
}

void
finalise_pooled_matrix (pooled_matrix pool, char_vector sptaxa)
{
  int i;
  finalise_spdist_matrix (pool->d_total);
  if (pool->d_total->n_missing) fprintf (stderr, "OBS: %d species pair combinations never appear on same gene familiy\n", pool->d_total->n_missing);
  for (i = 0; i < sptaxa->nstrings; i++) if (!pool->d_total->species_present[i])
    fprintf (stderr, "OBS: species \'%s\' never appears in data set; consider removing it from name list. Multifurcations also make this mess\n", sptaxa->string[i]);
  for (i = 0; i < pool->n_sptrees; i++) {
    finalise_spdist_matrix (pool->d[i]);
    complete_missing_spdist_from_global_spdist (pool->d[i], pool->d_total);
  }
}

void
find_maxtree_and_add_to_topol_space (spdist_matrix dist, pooled_matrix pool, topology_space tsp, int tree_method, bool use_within_gf_means)
{
  topology maxtree;
  copy_spdist_matrix_to_distance_matrix_upper (dist, pool->square_matrix, use_within_gf_means);
  maxtree = new_topology (tsp->taxlabel->nstrings);
  maxtree->taxlabel = tsp->taxlabel; tsp->taxlabel->ref_counter++; /* sptaxa is pointed here and at the calling function */
  if (tree_method == 0) bionj_from_distance_matrix (maxtree, pool->square_matrix); 
  else if (tree_method == 1) upgma_from_distance_matrix (maxtree, pool->square_matrix, false); // false -> upgma, true -> single linkage 
  else upgma_from_distance_matrix (maxtree, pool->square_matrix, true); // false -> upgma, true -> single linkage 
  add_topology_to_topology_space_if_distinct (maxtree, tsp, true); // true -> consider root location
  return;
}

