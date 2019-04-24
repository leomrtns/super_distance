#include "super_distance.h" 

void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .mode = arg_strn("m", "mode", "DNP", 0, 1, "distance (D), bipartition NJ (N), parsimony (P)"),
    .spname = arg_filen("s","species", "<species names>", 1, 1, "file with species names, one name per line (nexus-style bracketed comments are allowed)"),
    .outfil = arg_file0("o","output", "<newick>", "output file with species supertrees, in newick format (default '-')"),
    .genfil = arg_filen(NULL, NULL, NULL, 1, argc+2, "list of gene tree files, in newick format"),
    .end  = arg_end(20)
  };
  void* argtable[] = {params.help, params.mode, params.spname, params.outfil, params.genfil, params.end};
  params.argtable = argtable; 
  /* default values: */
  params.outfil->filename[0] = "-";
  params.mode->sval[0] = "DNP";
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
    printf ("Longer help explanation to come\n");
  }
  del_arg_parameters (params);
  exit (EXIT_FAILURE);
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
    sptrees = find_matrix_distance_species_tree (gene_nwk, species_names, 1e-7, false, false);
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
