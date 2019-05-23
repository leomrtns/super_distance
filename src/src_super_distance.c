#include "src_super_distance.h" 

void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .fast = arg_litn("F", "fast", 0, 1, "for too many leaves, estimates only two species trees"),
    .tol  = arg_dbl0("e", "epsilon", NULL, "tolerance (small value below which a branch length is considered zero for nodal distances)"),
    .spname = arg_file0("s","species", "<species names>", "file with species names, one per line (comments in square braces allowed); if absent, orthology is assumed"),
    .outfil = arg_file0("o","output", "<newick>", "output file with species supertrees, in newick format (default '-')"),
    .genfil = arg_filen(NULL, NULL, NULL, 1, argc+2, "list of gene tree files, in newick format"),
    .end  = arg_end(20)
  };
  void* argtable[] = {params.help, params.version, params.fast, params.tol, params.spname, params.outfil, params.genfil, params.end};
  params.argtable = argtable; 
  /* default values: */
  params.outfil->filename[0] = "-";
  params.tol->dval[0] = 1e-7;
  /* actual parsing: */
  if (arg_nullcheck(params.argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  if (arg_parse (argc, argv, params.argtable)) print_usage (params, argv[0]);
  //if (params.help->count) print_usage (params, argv[0]);
  return params;
}

void
del_arg_parameters (arg_parameters params)
{
  if (params.help) free (params.help);
  if (params.version) free (params.version);
  if (params.fast) free (params.fast);
  if (params.tol)  free (params.tol);
  if (params.spname) free (params.spname);
  if (params.outfil) free (params.outfil);
  if (params.genfil) free (params.genfil);
  if (params.end) free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.version->count) { printf ("%s\n", PACKAGE_VERSION); del_arg_parameters (params); exit (EXIT_SUCCESS); }

  if (params.end->count && (!params.help->count)) {  // params.end holds error messages
    arg_print_errors(stdout, params.end, basename(progname));
    printf ("Error when reading arguments from command line\n\n");
  }

  printf ("%s \n", PACKAGE_STRING);
  printf ("Matrix Representation with Distances: calculates pairwise distances between gene leaves, and estimates species trees from summary distance matrices\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("\nBased on several rescaled patristic distances, the program takes the average matrix between genes and estimates \n");
    printf (" the species tree using bioNJ, UPGMA and single-linkage after scaling back to the original values (more below). The program \n");
    printf (" also uses a distance matrix to project branch lengths on species trees missing lengths; \n\n");
    printf ("The branch length rescaling per gene can be the minimum, the average, the total sum, etc. and at the end these values\n");
    printf (" averaged over trees are scaled back in the final distance matrix, such that lengths in the supertree (species tree) are interpretable.\n");
    printf (" One exception is the nodal distance, which is based on the number of nodes between two leaves (e.g. NJst). In this case it may make\n");
    printf (" more sense to use another distance matrix to infer the branch lengths. Option 'F' uses averages distances projected on nodal-estimated tree; \n");
    printf (" it uses fewer scalings/options, providing a fast estimation. We avoid using individual gene trees since they may have \n");
    printf (" missing information (missing species or species pairs). For missing comparisons (when two species are never seen in the same gene tree)\n");
    printf (" we use the ultrametric condition (comparison to a common species) to estimate its value.\n\n");
    printf ("If a file with species names is given, the program allows for paralogs; otherwise it assumes orthology and that _at_least_ one tree has no missing data:\n");
    printf (" * Paralogy: the species names will be mapped to individual gene tree leaves (e.g. `ECOLI_a` and `ECOLI_b` will both map to species `ECOLI`).\n ");
    printf ("   Each gene tree can therefore have several copies of each species, and can also have missing species.\n");
    printf (" * Orthology: if a file with species names is not given, however, it is assumed that each species is represented at most once per gene, and\n");
    printf ("   furthermore that the leaf names represent the species, and are thus identical across trees. This mode is the underlying assumption behind\n");
    printf ("   most tree comparison software, although here missing data for some trees (not all) is allowed. I.e. as long as one tree has full information\n");
    printf ("   (for all species), then others can have some absent species.\n");
    printf ("With paralogs or not, it is not recommended to have missing entries in the distance matrix (e.g. when a species pair does not appear in any tree),\n");
    printf (" and matrix representation with distances methods work better with more 'complete' gene trees.\n\n");
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i;
  char_vector original_spnames = NULL, species_names = NULL;
  newick_space gene_nwk = new_newick_space();
  newick_space sptrees = NULL;
  FILE *stream=NULL;
  char *s;

  /* initialisation: */
  biomcmc_random_number_init (0);
  arg_parameters params = get_parameters_from_argv (argc, argv);

  for (i=0; i < params.genfil->count; i++) update_newick_space_from_file (gene_nwk, (char*) params.genfil->filename[i]);
  /* read species names and reduce to valid ones: */ 
  if (params.spname->count) {
    original_spnames = new_char_vector_from_file ((char*) params.spname->filename[0]);
    species_names = get_species_names_from_newick_space (gene_nwk, original_spnames, true); // true= reorder species names
    del_char_vector (original_spnames); // free() or decrease ref_counter is same as species_names
  }
  else 
    species_names = assume_species_names_from_newick_space (gene_nwk); // true= reorder species names

  /* open output file: */
  if (strstr (params.outfil->filename[0], "-")) stream = stdout;
  else stream = biomcmc_fopen (params.outfil->filename[0], "w");

  /* mode: fast version of patristic-based upgma tree: */
  if (params.fast->count) { 
    sptrees = find_matrix_distance_species_tree (gene_nwk, species_names, params.tol->dval[0], false, false, true);
    for (i=0; i < sptrees->ntrees; i++) { 
      s = topology_to_string_by_name (sptrees->t[i], sptrees->t[i]->blength);
      fprintf (stream, "[F%02d] %s\n", i, s); fflush(stream); free (s);
    }
    del_newick_space (sptrees);
  }
  /* mode: patristic-based nj/upgma tree: */
  else { 
    sptrees = find_matrix_distance_species_tree (gene_nwk, species_names, params.tol->dval[0], false, false, false);
    for (i=0; i < sptrees->ntrees; i++) { 
      s = topology_to_string_by_name (sptrees->t[i], sptrees->t[i]->blength);
      fprintf (stream, "[D%02d] %s\n", i, s); fflush(stream); free (s);
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
