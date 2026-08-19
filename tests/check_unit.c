#include <super_distance.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

int
main ()
{
  int i, j, failed = 0;
  char ortholog_tree1[] = "(((Alpha:1,Beta:2):1,(Gamma:2,Delta:1):1):2,((Epsilon:1,Zeta:2):1,Eta:2):1)";
  char ortholog_tree2[] = "(((Alpha:2,Beta:1):1,(Gamma:1,Delta:2):1):1,(Epsilon:2,Theta:1):2)";
  char ortholog_tree3[] = "(Zeta:2,Eta:1,Theta:2)";
  char paralog_tree[] = "((((Alpha_1:1,Alpha_2:2):1,(Beta:1,Gamma:2):1):1,((Delta:1,Epsilon:2):1,(Zeta:2,Eta:1):1):1):1,Theta:2)";
  const char *expected[] = {"Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta"};
  bool found[8] = {false, false, false, false, false, false, false, false};
  newick_space orthologs = new_newick_space (), paralogs = new_newick_space (), sptrees;
  char_vector species, supplied_species = new_char_vector (8), paralog_species;

  biomcmc_random_number_init (0);
  update_newick_space_from_string (orthologs, ortholog_tree1, strlen (ortholog_tree1));
  update_newick_space_from_string (orthologs, ortholog_tree2, strlen (ortholog_tree2));
  update_newick_space_from_string (orthologs, ortholog_tree3, strlen (ortholog_tree3));
  species = assume_species_names_from_newick_space (orthologs);
  for (i = 0; i < species->nstrings; i++) for (j = 0; j < 8; j++)
    if (!strcmp (species->string[i], expected[j])) found[j] = true;
  if (species->nstrings != 8) failed = 1;
  for (i = 0; i < 8; i++) if (!found[i]) failed = 1;

  sptrees = find_matrix_distance_species_tree (orthologs, species, 1e-7, false, false, false, true);
  if (sptrees->ntrees != 18) failed = 1;
  del_newick_space (sptrees);

  sptrees = find_matrix_distance_species_tree (orthologs, species, 1e-7, false, false, false, false);
  if (sptrees->ntrees != 18) failed = 1;
  del_newick_space (sptrees);

  for (i = 0; i < 8; i++) {
    char_vector_add_string (supplied_species, expected[i]);
  }
  update_newick_space_from_string (paralogs, paralog_tree, strlen (paralog_tree));
  paralog_species = get_species_names_from_newick_space (paralogs, supplied_species, true);
  sptrees = find_matrix_distance_species_tree (paralogs, paralog_species, 1e-7, false, false, false, false);
  if (sptrees->ntrees != 36) failed = 1;

  del_newick_space (sptrees);
  del_char_vector (paralog_species);
  del_char_vector (supplied_species);
  del_newick_space (paralogs);
  del_char_vector (species);
  del_newick_space (orthologs);
  biomcmc_random_number_finalize ();
  return failed ? TEST_FAILURE : TEST_SUCCESS;
}
