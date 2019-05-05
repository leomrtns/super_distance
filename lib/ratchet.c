typedef struct genetree_struct* genetree;
typedef struct gene_sptrees_struct* gene_sptrees;

struct genetree_struct
{
  topology tree;
  double minmax[2 * NDISTS], dist[NDISTS]; // values are integers; might have a variable "scale" to allow for quick or no scaling 
  splitset split;
};

struct gene_sptrees_struct
{
  int n_genes, n_ratchet, n_proposal;
  int next_avail; // idx to ratchet elements 
  bool optimise_locally;
  genetree *gene;
  topology *ratchet, *proposal; // may need vector of distances
  double *ratchet_score, best_score; // I'm forgetting other vars like minmax etc. but for now we will use simple score (gene dists are rescaled)
};


int
main (int argc, char **argv)
{

  /* k samplings of i-percent above of genefams and search for optimal trees, using previous sptrees if exist */
  for (k = 0; k < n_samples; k++) {
    /* chose a new set of genefams */ 
    initialise_gene_sptrees_with_genetree_filenames (gs, genefiles);
    /* few rounds of optimisation using global and local scores */
    for (j = 0; j < n_cycles; j++) { 
      /* optimisation of ratchet trees */
      if (!j) {
        sorting_of_gene_sptrees_ratchet (gs, tsp, true); // bool decides if local score or not, can be (bool)j%2
        improve_gene_sptrees_ratchet (gs, 1); // spend just one iteration with local score trees
      }
      else {
        sorting_of_gene_sptrees_ratchet (gs, tsp, false); // bool decides if local score or not, can be (bool)j%2
        improve_gene_sptrees_ratchet (gs, n_iter_per_cycle);// multiples of 4 are better since there are 4 randomisers
      }
      if (tsp) { del_topology_space (tsp); tsp = NULL; } // discard after first use
    }
    /* save top best trees to file */
    for (i = 0; i < n_output_trees; i++) { // i must be smaller than gs->n_ratchet
      j = (i + gs->next_avail+1) % gs->n_ratchet; // circular index 
      s = topology_to_string_by_name (gs->ratchet[j], NULL); 
      fprintf (stream, "tree PAUP_%d = %s;\n", k * n_output_trees + i, s); fflush(stream); free (s); // each k sampling will have i x k trees
    }
    /* print best tree to screen */
    j = (gs->next_avail+1) % gs->n_ratchet; // circular index 
    s = topology_to_string_by_name (gs->ratchet[j], NULL); /* second parameter is vector with branch lengths */
    printf ("[sampling %d score %lf] %s;\n",k, gs->ratchet_score[j],s); fflush(stdout); free (s);
    time1 = clock (); fprintf (stderr, "partial timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC); time0 = time1;
  } // for (k<number of subsamples)
  fprintf (stream, "End;\n");
  fclose (stream);

	return EXIT_SUCCESS;
}

/* genetree_struct functions */

genetree 
new_genetree (topology g_tree, topology s_tree)
{
  int i;
  genetree gt;
  gt = (genetree) biomcmc_malloc (sizeof (struct genetree_struct));
  gt->tree = g_tree; gt->tree->ref_counter++; // point to element of topology_space (usually)

  gt->split = create_splitset_dSPR_genespecies (gt->tree, s_tree);
  // init_tree_recon_from_species_topology() would order sptree names everytime (and already calculates distances...)
  if (!gt->tree->rec) gt->tree->rec = new_reconciliation (g_tree->nleaves, s_tree->nleaves); // if gene has rec, it better be from same sptree size!
  index_sptaxa_to_reconciliation (s_tree->taxlabel, g_tree->taxlabel, gt->tree->rec);
  for (i = 0; i < NDISTS; i++) {
    gt->minmax[i] = 1.e35;  // higher bound for min
    gt->minmax[i+NDISTS] = -1.e35; // lower bound for max
    gt->dist[i] = -1.;
  }
//  del_char_vector (g_tree->taxlabel); // should just decrease ref_counter and then delete later, with topology_space
  return gt;
}

void
del_genetree (genetree gt)
{
  if (!gt) return;
  del_splitset (gt->split);
  del_topology (gt->tree);
  free (gt);
  return;
}

double
calculate_genetree_distance (genetree gt, topology s_tree, bool local_optimum) 
{ /* local_optima -> return zero if a new local optimum is found (new minimum for _some_ distance)  */
  int k;
  double g_score = 0.;
  gene_tree_reconcile (gt->tree, s_tree);
//  dSPR_gene_species_rf (gt->tree, s_tree, gt->split); // does not compute hdist and dSPR, only RF
  dSPR_gene_species (gt->tree, s_tree, gt->split); 
  gt->dist[0] = (double) gt->tree->rec->ndups;
  gt->dist[1] = (double) gt->tree->rec->nloss;
  gt->dist[2] = (double) gt->tree->rec->ndcos;
  gt->dist[3] = (double) gt->split->rf;
  gt->dist[4] = (double) gt->split->hdist;
  gt->dist[5] = (double) gt->split->spr + gt->split->spr_extra;  // order same as guenomu, NOT same as genefam_dist or treesignal.py
  /* minimum values are always updated */
  for (k=0; k < NDISTS; k++) if (gt->minmax[k] > gt->dist[k]) gt->minmax[k] = gt->dist[k]; // min
  if (local_optimum) {
    g_score = 1.;  // score is 0 if at least one min is found, 1.1 if one new max is found, and 1 o.w.
    for (k=0; k < NDISTS; k++) if (gt->minmax[k] == gt->dist[k]) g_score = 0.;
    for (k=0; k < NDISTS; k++) { if (gt->minmax[k+NDISTS] < gt->dist[k]) { gt->minmax[k+NDISTS] = gt->dist[k]; g_score = 1.1; }} // max
  }
  else for (k=0; k < NDISTS; k++) g_score += biomcmc_log1p( (double)(gt->dist[k])/(double)(gt->minmax[k+NDISTS])); 
  // else  for (k=0; k < NDISTS; k++) g_score += biomcmc_log1p( (double)(gt->dist[k])); 
  // for (k=0; k < NDISTS; k++) printf ("%.0f (%.0f) %.0f | ", gt->minmax[k], gt->dist[k], gt->minmax[k+NDISTS]); // DEBUG
  return g_score;
}

gene_sptrees
new_gene_sptrees (int n_genes, int n_ratchet, int n_proposal, topology_space sptree)
{
  int i;
  gene_sptrees gs;
  gs = (gene_sptrees) biomcmc_malloc (sizeof (struct gene_sptrees_struct));
  gs->n_genes = n_genes;
  gs->n_ratchet = n_ratchet;
  gs->n_proposal = n_proposal;
  gs->optimise_locally = false;
  gs->best_score = 1.e35;
  if (gs->n_ratchet < 10) gs->n_ratchet = 10; // ratchet[0] has best score at first (due to initial_sorting() )
  if (gs->n_proposal < 2) gs->n_proposal = 2;
  gs->next_avail = gs->n_ratchet - 1; // idx of currently worse tree (which is just before best tree in a ratchet)

  gs->gene = (genetree*) biomcmc_malloc (gs->n_genes * sizeof (genetree)); 
  gs->ratchet  = (topology*) biomcmc_malloc (gs->n_ratchet * sizeof (topology)); 
  gs->proposal = (topology*) biomcmc_malloc (gs->n_proposal * sizeof (topology));
  gs->ratchet_score  = (double*) biomcmc_malloc (gs->n_ratchet * sizeof (double)); 

  for (i=0; i < gs->n_ratchet; i++) {
    gs->ratchet[i] = init_topol_for_new_gene_sptrees (sptree->distinct[0]);
    randomize_topology (gs->ratchet[i]); // should be filled outside this function (to recycle optimal sptrees from previous iteration)
  }
  for (i=0; i < gs->n_proposal; i++) gs->proposal[i] = init_topol_for_new_gene_sptrees (sptree->distinct[0]);
  for (i=0; i < gs->n_genes; i++) gs->gene[i] = NULL; // must be created later, when we receive the gene topol_spaces 

  initialise_ratchet_from_topol_space (gs, sptree);
  return gs;
}

topology
init_topol_for_new_gene_sptrees (topology original)
{
  topology tree = new_topology(original->nleaves);
  tree->taxlabel = original->taxlabel;
  tree->taxlabel->ref_counter++;
  new_mrca_for_topology (tree);
  return tree;
}

void
initialise_ratchet_from_topol_space (gene_sptrees gs, topology_space sptree)
{
  int i, n_trees_to_copy = sptree->ndistinct;
  if (n_trees_to_copy > gs->n_ratchet) n_trees_to_copy = gs->n_ratchet;
  for (i=0; i < n_trees_to_copy; i++) copy_topology_from_topology (gs->ratchet[gs->n_ratchet - i -1], sptree->distinct[i]); 
  // the remaining trees (from n_trees_to_copy to n_ratchet) are random (if it's first time this function is called) 
  // or are best from previous iteration (first on ratchet) 
  return;
}

void
initialise_gene_sptrees_with_genetree_filenames (gene_sptrees gs, char_vector gene_files)
{ // sample gene families, and calculate max values 
  int i, j, *idx;
  topology_space genetre;
  genetree pivot;

  idx = (int*) biomcmc_malloc (gene_files->nstrings * sizeof (int)); 
  for (i=0; i < gene_files->nstrings; i++) idx[i] = i;

  for (i=0; i < gs->n_genes; i++) {
    j = biomcmc_rng_unif_int (gene_files->nstrings - i);
    genetre = read_topology_space_from_file (gene_files->string[ idx[j] ], NULL, false);// read i-th gene family; false means that rooting is neglected 
    idx[j] = idx[gene_files->nstrings - i -1]; // avoids replacement
    pivot = gs->gene[i];
    gs->gene[i] = new_genetree (genetre->distinct[0], gs->proposal[0]);
    for (j = 0; j < 32; j++) {
      randomize_topology (gs->proposal[0]);
      calculate_genetree_distance (gs->gene[i], gs->proposal[0], true); // true --> local (since we must find max values)
    }
    del_genetree (pivot);
    del_topology_space (genetre);
  }
  if (idx) free (idx);
  return;
}

void
del_gene_sptrees (gene_sptrees gs)
{
  int i;
  if (!gs) return;
  if (gs->gene) {
    for (i = gs->n_genes - 1; i >= 0 ; i--) del_genetree (gs->gene[i]);
    free (gs->gene);
  }
  if (gs->ratchet) {
    for (i = gs->n_ratchet - 1; i >= 0 ; i--) del_topology (gs->ratchet[i]);
    free (gs->ratchet);
  }
  if (gs->proposal) {
    for (i = gs->n_proposal - 1; i >= 0 ; i--) del_topology (gs->proposal[i]);
    free (gs->proposal);
  }
  if (gs->ratchet_score) free (gs->ratchet_score);
  free (gs);
  return;
}

void
sorting_of_gene_sptrees_ratchet (gene_sptrees gs, topology_space tsp, bool local_optimum)
{ // first step before start optimisation: sort ratchet 
  topology *pivot; // temporary vector with original location of topologies
  empfreq_double ef;
  double *score_vec = NULL;
  gs->optimise_locally = local_optimum;
  int i, min_size;

  /* make sure ratchet has best scores (local scores initially, to update MinMax), but not necessarily in order */
  if (tsp) {
    min_size = (tsp->ndistinct < gs->n_ratchet)? tsp->ndistinct : gs->n_ratchet;
    score_vec  = (double*) biomcmc_malloc (tsp->ndistinct * sizeof (double)); 
    for (i=0; i < tsp->ndistinct; i++) score_vec[i] = calculate_species_tree_score (gs, tsp->distinct[i], true);
    /* Find the score of each tree and sort their scores from lowest (best to highest (worst)  */ 
    ef = new_empfreq_double_sort_increasing (score_vec, tsp->ndistinct);
    for (i=0; i < min_size; i++) copy_topology_from_topology (gs->ratchet[i], tsp->distinct[ ef->d[i].idx ]);
    for (i=0; i < min_size; i++) gs->ratchet_score[i] = ef->d[i].freq; // doesnt need pivot since ef->freqs came from ratchet_score[]
    if (score_vec) free (score_vec);
  }
  else min_size = 0; 
  /* first (or all) scores come from topol_space, remaining (or all) come from previous ratchet/sample (that are initialised to random trees) */
  for (i = min_size; i < gs->n_ratchet; i++) gs->ratchet_score[i] = calculate_species_tree_score (gs, gs->ratchet[i], true);
  /* working with local scores so far, recalculation may be needed */
  if (!gs->optimise_locally) for (i=0; i < gs->n_ratchet; i++) gs->ratchet_score[i] = calculate_species_tree_score (gs, gs->ratchet[i], gs->optimise_locally);

  improve_gene_sptrees_initial_state (gs); // try to replace very bad initial trees by better ones

  ef = new_empfreq_double_sort_increasing (gs->ratchet_score, gs->n_ratchet);
  pivot = (topology*) biomcmc_malloc (gs->n_ratchet * sizeof (topology)); 
  for (i=0; i < gs->n_ratchet; i++) pivot[i] = gs->ratchet[i]; //vector version of tmp=a;a=b;b=tmp
  for (i=0; i < gs->n_ratchet; i++) gs->ratchet[i] = pivot[ ef->d[i].idx ]; 
  for (i=0; i < gs->n_ratchet; i++) gs->ratchet_score[i] = ef->d[i].freq; // doesnt need pivot since ef->freqs came from ratchet_score[]
  gs->best_score = ef->d[0].freq; // lowest (best) score is first element (may need to be recalc every time new proposals are created, using new minima)
  fprintf (stderr, "Ordered list of best scores: ");
  for (i = 0; i < 5; i++) printf ("%7.3lf ", gs->ratchet_score[i]);   printf (" ...... ");
  for (i = gs->n_ratchet - 6; i < gs->n_ratchet; i++) printf ("%7.3lf ", gs->ratchet_score[i]);   printf ("\n");
  gs->next_avail = gs->n_ratchet - 1; // idx of currently worse tree (which is just before best tree in a ratchet)

  del_empfreq_double (ef);
  if (pivot) free (pivot);
  return;
}

void // OLD
sorting_of_gene_sptrees_ratchet_without_optimisation (gene_sptrees gs, topology_space tsp, bool local_optimum)
{ // first step before start optimisation: sort ratchet 
  topology *pivot; // temporary vector with original location of topologies
  empfreq_double ef;
  double *score_vec = NULL;
  gs->optimise_locally = local_optimum;
  int i;

  if ((!tsp) || (tsp->ndistinct <= gs->n_ratchet)) {
    /* calculate twice to ensure MinMax are updated */
    for (i=0; i < gs->n_ratchet; i++) gs->ratchet_score[i] = calculate_species_tree_score (gs, gs->ratchet[i], true);
    for (i=0; i < gs->n_ratchet; i++) gs->ratchet_score[i] = calculate_species_tree_score (gs, gs->ratchet[i], gs->optimise_locally);
    /* Find the score of each tree and sort their scores from lowest (best to highest (worst)  */ 
    ef = new_empfreq_double_sort_increasing (gs->ratchet_score, gs->n_ratchet);
    pivot = (topology*) biomcmc_malloc (gs->n_ratchet * sizeof (topology)); 
    /* reorder ratchet, by creating temporary pointers with previous location */
    for (i=0; i < gs->n_ratchet; i++) pivot[i] = gs->ratchet[i]; //vector version of tmp=a;a=b;b=tmp
    for (i=0; i < gs->n_ratchet; i++) gs->ratchet[i] = pivot[ ef->d[i].idx ]; 
    if (pivot) free (pivot);
  }
  else {
    score_vec  = (double*) biomcmc_malloc (tsp->ndistinct * sizeof (double)); 
    /* calculate twice to ensure MinMax are updated */
    for (i=0; i < tsp->ndistinct; i++) score_vec[i] = calculate_species_tree_score (gs, tsp->distinct[i], true);
    for (i=0; i < tsp->ndistinct; i++) score_vec[i] = calculate_species_tree_score (gs, tsp->distinct[i], gs->optimise_locally);
    /* Find the score of each tree and sort their scores from lowest (best to highest (worst)  */ 
    ef = new_empfreq_double_sort_increasing (score_vec, tsp->ndistinct);
    for (i=0; i < gs->n_ratchet; i++) copy_topology_from_topology (gs->ratchet[i], tsp->distinct[ ef->d[i].idx ]);
    if (score_vec) free (score_vec);
  }

  for (i=0; i < gs->n_ratchet; i++) gs->ratchet_score[i] = ef->d[i].freq; // doesnt need pivot since ef->freqs came from ratchet_score[]
  gs->best_score = ef->d[0].freq; // lowest (best) score is first element (may need to be recalc every time new proposals are created, using new minima)
  fprintf (stderr, "Ordered list of best scores: ");
  for (i = 0; i < 5; i++) printf ("%7.3lf ", gs->ratchet_score[i]);   printf (" ...... ");
  for (i = gs->n_ratchet - 6; i < gs->n_ratchet; i++) printf ("%7.3lf ", gs->ratchet_score[i]);   printf ("\n");
  gs->next_avail = gs->n_ratchet - 1; // idx of currently worse tree (which is just before best tree in a ratchet)

  del_empfreq_double (ef);
  return;
}

double
calculate_species_tree_score (gene_sptrees gs, topology s_tree, bool local_optimum)
{
  int i;
  double score = 0.;
  for (i = 0; i < gs->n_genes; i++) score += calculate_genetree_distance (gs->gene[i], s_tree, local_optimum);
  return score;
} 

void
improve_gene_sptrees_initial_state (gene_sptrees gs)
{ 
  int i, j, iteration = 0; // iteration will determine the sets of possible branch swaps 
  double score = 0;
  for (j = 0; j < gs->n_ratchet; j++) { 
    generate_new_gene_sptrees_proposal_trees (gs, j, iteration); // this function guarantees that proposal trees not the same as any in ratchet 
    for (i = 0; i < gs->n_proposal; i++) {
      score = calculate_species_tree_score (gs, gs->proposal[i], gs->optimise_locally);
      if (score < gs->ratchet_score[j]) {
        topology pivot = gs->proposal[i];
        gs->proposal[i] = gs->ratchet[j];
        gs->ratchet[j] = pivot;
        gs->ratchet_score[j] = score;
      }
    }
  }
  return;
}

void 
improve_gene_sptrees_ratchet (gene_sptrees gs, int n_iterations)
{ 
  int i,j, iteration;
  double score = 0;
  if (n_iterations < 1) n_iterations = 1;
  fprintf (stderr, "Best score at end of each iteration: ");
  for (iteration = 0; iteration < n_iterations; iteration++) {
    for (j = 0; j < gs->n_ratchet; j++) { 
      generate_new_gene_sptrees_proposal_trees (gs, j, iteration);
      for (i = 0; i < gs->n_proposal; i++) {
        score = calculate_species_tree_score (gs, gs->proposal[i], gs->optimise_locally);
        if (score <= gs->best_score) add_topol_to_gene_sptrees_ratchet (gs, i, score);
      }
    }
    fprintf (stderr, "%7.3lf ", gs->ratchet_score[ (gs->next_avail+1) % gs->n_ratchet ]);
  } // for (iteration)
  fprintf (stderr, "\n");
  return;
}

void
generate_new_gene_sptrees_proposal_trees (gene_sptrees gs, int idx, int iteration)
{
  int coinflip = 0;
  int i, j;
  for (i = 0; i < gs->n_proposal; i++) copy_topology_from_topology (gs->proposal[i], gs->ratchet[idx]);
  if (!(iteration%2)) {
    for (i = 0; i < gs->n_proposal; i++) {
      topology_apply_rerooting (gs->proposal[i], false);
      for (j = 0; j < biomcmc_rng_unif_int (3); j++) topology_apply_shortspr (gs->proposal[i], false);
      for (j = 0; j < biomcmc_rng_unif_int (2); j++) topology_apply_nni (gs->proposal[i], false);
    }
  }
  if (!((iteration+1)%2)) {
    for (i = 0; i < gs->n_proposal; i++) {
      for (j = 0; j <= biomcmc_rng_unif_int (3); j++) topology_apply_nni (gs->proposal[i], false);
      for (j = 0; j <  biomcmc_rng_unif_int (3); j++) topology_apply_spr_unrooted (gs->proposal[i], false);
    }
  }
  // if rooted trees are the same, randomise one of them (assuming random tree won't be the same anymore...)
  for (i = 0; i < gs->n_proposal; i++) {
    for (j = 0; j < i; j++) if (topology_is_equal (gs->proposal[i], gs->proposal[j])) { randomize_topology (gs->proposal[i]); j = i; }
    if (j < i) for (j = 0; j < gs->n_ratchet; j++) if (topology_is_equal (gs->proposal[i], gs->ratchet[j])) { randomize_topology (gs->proposal[i]); j = gs->n_ratchet; }
  }
  return;
}

void // OLD
generate_new_gene_sptrees_proposal_trees_deterministic (gene_sptrees gs, int idx, int iteration)
{
  int coinflip = 0;
  int i, j, n_2 = gs->n_proposal/2;
  for (i = 0; i < gs->n_proposal; i++) copy_topology_from_topology (gs->proposal[i], gs->ratchet[idx]);
  if (!(iteration%4)) {
    for (i = 0; i < n_2; i++) topology_apply_rerooting (gs->proposal[i], false);
    for (i = n_2; i < gs->n_proposal; i++) topology_apply_shortspr (gs->proposal[i], false);
  }
  if (!((iteration+1)%4)) {
    for (i = 0; i < n_2; i++) topology_apply_spr_unrooted (gs->proposal[i], false);
    for (i = n_2; i < gs->n_proposal; i++) topology_apply_nni (gs->proposal[i], false);
  }
  if (!((iteration+2)%4)) {
    for (i = 0; i < gs->n_proposal; i++) topology_apply_spr_unrooted (gs->proposal[i], false);
    for (i = 0; i < n_2; i++) topology_apply_nni (gs->proposal[i], false);
  }
  if (!((iteration+3)%4)) {
    for (i = 0; i < gs->n_proposal; i++) { 
      topology_apply_shortspr (gs->proposal[i], false); 
      topology_apply_spr_unrooted (gs->proposal[i], false); 
    }
  }
  // if rooted trees are the same, randomise one of them (assuming random tree won't be the same anymore...)
  for (i = 0; i < gs->n_proposal; i++) for (j = 0; j < i; j++) if (topology_is_equal (gs->proposal[i], gs->proposal[j])) randomize_topology (gs->proposal[j]); 
  for (i = 0; i < gs->n_proposal; i++) {
    coinflip = biomcmc_rng_unif_int (10); 
    if (coinflip < 2) topology_apply_spr_unrooted (gs->proposal[i], false); // 2/10 chance
    if (coinflip > 8) randomize_topology (gs->proposal[i]); // 1/10 chance
  }
  return;
}

void
add_topol_to_gene_sptrees_ratchet (gene_sptrees gs, int idx_proposal, double score)
{
  topology pivot = gs->proposal[idx_proposal];
  gs->proposal[idx_proposal] = gs->ratchet[gs->next_avail];
  gs->ratchet[gs->next_avail] = pivot;
  gs->ratchet_score[gs->next_avail] = score;
  // printf ("DEBUG::bestfound::%d  \t%lf\n",gs->next_avail, score);
  gs->next_avail--;
  if (gs->next_avail < 0) gs->next_avail = gs->n_ratchet - 1;
  gs->best_score = score;
  return;
}
