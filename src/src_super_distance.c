#include "src_super_distance.h" 

void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .mode = arg_strn("m", "mode", "DNP", 0, 1, "distance (D), bipartition NJ (N), parsimony (P) <only D currently functional>"),
    .tol  = arg_dbl0("e", "epsilon", NULL, "tolerance (small value below which a branch length is considered zero for nodal distances)"),
    .spname = arg_filen("s","species", "<species names>", 1, 1, "file with species names, one name per line (nexus-style bracketed comments are allowed)"),
    .outfil = arg_file0("o","output", "<newick>", "output file with species supertrees, in newick format (default '-')"),
    .genfil = arg_filen(NULL, NULL, NULL, 1, argc+2, "list of gene tree files, in newick format"),
    .end  = arg_end(20)
  };
  void* argtable[] = {params.help, params.mode, params.tol, params.spname, params.outfil, params.genfil, params.end};
  params.argtable = argtable; 
  /* default values: */
  params.outfil->filename[0] = "-";
  params.mode->sval[0] = "D";
  params.tol->dval[0] = 1e-7;
  /* actual parsing: */
  if (arg_nullcheck(params.argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  if (arg_parse (argc, argv, params.argtable)) print_usage (params, argv[0]);
  if (params.help->count) print_usage (params, argv[0]);
  return params;
}

void
del_arg_parameters (arg_parameters params)
{
  if (params.help) free (params.help);
  if (params.mode) free (params.mode);
  if (params.tol)  free (params.tol);
  if (params.spname) free (params.spname);
  if (params.outfil) free (params.outfil);
  if (params.genfil) free (params.genfil);
  if (params.end) free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.end->count && (!params.help->count)) {  // params.end holds error messages
    arg_print_errors(stdout, params.end, basename(progname));
    printf ("Error when reading arguments from command line\n\n");
  }
  printf ("Calculates pairwise distances between gene leaves, and estimates species trees from summary distance matrices\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("Currently only the distance-baseed approaches are implemented:\n");
    printf ("Based on several rescaled patristic distances, the program takes the average matrix between genes and estimates \n");
    printf (" the species tree using bioNJ, UPGMA and single-linkage after scaling back to the original values (more below). The program \n");
    printf (" also can use an 'external' distance matrix and project branch lengths on it; we are now playing with these options \n");
    printf (" and soon the user will have more control over which to output\n\n");
    printf ("The branch length rescaling per gene can be the minimum, the average, the total sum, etc. and at the end these values\n");
    printf (" averaged over trees are scaled back in the final distance matrix, such that the supertree (species tree) lengths are interpretable.\n");
    printf (" One exception is the nodal distance, which is based on the number of nodes between two leaves (e.g. NJst). In this case it may make\n");
    printf (" more sense to use another distance matrix to infer the branch lengths. We avoid using individual gene trees since they may have \n");
    printf (" missing information (missing species or species pairs). For missing comparisons (when two species are never seen in the same gene tree)\n");
    printf (" we use the ultrametric condition (comparison to a common species) to estimate its value.\n");
  }
  del_arg_parameters (params);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i;
  char_vector original_spnames, species_names;
  newick_space gene_nwk = new_newick_space();
  newick_space sptrees;
  FILE *stream=NULL;
  char *s;

  /* initialisation: */
  biomcmc_random_number_init (0);
  arg_parameters params = get_parameters_from_argv (argc, argv);

  /* read species names and reduce to valid ones: */ 
  original_spnames = new_char_vector_from_file ((char*) params.spname->filename[0]);
  for (i=0; i < params.genfil->count; i++) update_newick_space_from_file (gene_nwk, (char*) params.genfil->filename[i]);
  species_names = get_species_names_from_newick_space (gene_nwk, original_spnames, true); // true= reorder species names
  del_char_vector (original_spnames); // free() or decrease ref_counter is same as species_names

  /* open output file: */
  if (strstr (params.outfil->filename[0], "-")) stream = stdout;
  else stream = biomcmc_fopen (params.outfil->filename[0], "w");

  /* mode: patristic-based nj/upgma tree: */
  if (strstr(params.mode->sval[0],"D")) {
    sptrees = find_matrix_distance_species_tree (gene_nwk, species_names, params.tol->dval[0], false, false);
    for (i=0; i < sptrees->ntrees; i++) { 
      s = topology_to_string_by_name (sptrees->t[i], sptrees->t[i]->blength);
      fprintf (stream, "[D%02d] %s\n", i, s); fflush(stream); free (s);
    }
    del_newick_space (sptrees);
  }
  /* mode: nj/upgma tree based on bipartition (MRP) distances; usually initial MRP tree: */
  if (strstr(params.mode->sval[0],"N")) {
    sptrees = find_upgma_mrp_species_tree (gene_nwk, species_names, false, false);
    for (i=0; i < sptrees->ntrees; i++) { 
      s = topology_to_string_by_name (sptrees->t[i], sptrees->t[i]->blength);
      fprintf (stream, "[N%02d] %s\n", i, s); fflush(stream); free (s);
    }
    del_newick_space (sptrees);
  }

  if (stream != stdout) fclose (stream);

  biomcmc_random_number_finalize ();
  del_newick_space (gene_nwk);
  del_newick_space (sptrees);
  del_char_vector (species_names);
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}
