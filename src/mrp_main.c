#include <super_sptree.h> 


void print_welcome (char *progname);
void print_usage (char *progname);

void
print_welcome (char *progname)
{
  fprintf (stderr, "  Welcome to %s -- Estimation of supertree from gene trees (*not* multi-labelled) using Matrix Representation with Parsimony\n\n", basename (progname));
}

void
print_usage (char *progname)
{
  fprintf (stderr, "\tUSAGE: %s <species file> <gene tree 1> ... <gene tree N>\n\n", basename (progname));
  fprintf (stderr, "where  < species file >      is the name of a text file with the names of the species\n");
  fprintf (stderr, "       < gene tree ?? >      are file names for the gene family trees\n");
  fprintf (stderr, "\n  - remember that the species names should be found within the gene names\n");
  exit (EXIT_FAILURE);
}

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  char_vector species;
  char *s;
  int i, score, *idx_gene_to_sp = NULL;
  topology_space genetre;
  topology sptree;
  binary_parsimony pars;
  empfreq ef;

  time0 = clock ();
  print_welcome (argv[0]);
  if (argc < 3) print_usage (argv[0]);
  biomcmc_random_number_init(0);

  species = new_char_vector_from_file (argv[1]);
  char_vector_remove_duplicate_strings (species); /* duplicate names */
  char_vector_longer_first_order (species); // order species names from longest to shortest (speed up gene/spnames comparison) 
  ef = new_empfreq_sort_decreasing (species->nchars, species->nstrings, 1); /* 1=size_t (0=char, 2=int)*/

  pars = new_binary_parsimony (species->nstrings); /* no sptree information, just a matrix from genetrees */

  for (i = 2; i < argc; i++) { /* scan through gene files */
    genetre  = read_topology_space_from_file (argv[i], NULL, false);// read i-th gene family; false -> neglects rooting
    idx_gene_to_sp = (int *) biomcmc_realloc ((int*) idx_gene_to_sp, genetre->distinct[0]->nleaves * sizeof (int));/* for each gene leaf, index of species */
    index_sptaxa_to_genetaxa (species, genetre->taxlabel, idx_gene_to_sp, ef);/* map species names to gene names and store into idx[] */
    update_binary_parsimony_from_topology (pars, genetre->distinct[0], idx_gene_to_sp);
    del_topology_space (genetre); 
  }

  sptree = new_topology (species->nstrings);
  sptree->taxlabel = species; sptree->taxlabel->ref_counter++; 
  randomise_topology (sptree);
  new_mrca_for_topology (sptree);
  
  time1 = clock (); fprintf (stderr, "timing: %.8f secs\n\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC); time0 = time1;

  for (i=0; i < 10; i++) {
    topology_apply_spr_unrooted (sptree, false);
    score = binary_mrp_parsimony_score_of_topology (mrp, sptree);
    s = topology_to_string_by_name (sptree, NULL);
    printf (" %6d\t %s; \n", score, s); 
  }

  del_char_vector (species);
  if (idx_gene_to_sp) free (idx_gene_to_sp);
  del_empfreq (ef);
  del_mrp_parsimony (mrp);
  del_topology (sptree);
  biomcmc_random_number_finalize();

	time1 = clock (); fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

	return EXIT_SUCCESS;
}

