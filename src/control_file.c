#include "control_file.h" 

control_file
new_control_file_from_argv (int argc, char **argv)
{
  int i;
  struct arg_lit  *at_help = arg_lit0("h","help","print a longer help and exit");
  struct arg_dbl  *at_lbda = arg_dbln("L","lambda", NULL, 0, Ndist, "baseline reconciliation hyperprior that control distances (one or several values, should be << 1)");
  struct arg_int  *at_nsa1 = arg_int0("A","sa_samples",NULL, "number of cycles/samples simulated Annealing");
  struct arg_str  *at_dstr = arg_str0("D","distances",NULL, "string with seven  1's and 0's describing distances to be used [dup, los, ils, rf, hdist, hdist_final, spr]");
  struct arg_file *at_spec = arg_file0("S","species" ,NULL, "file with species names, one name per line (nexus-style bracketed comments are allowed)");
  struct arg_end  *at_end   = arg_end(20);
  void* argtable[] = {at_help, at_lbda, at_nsa1, at_dstr, at_spec, at_end};

  if (arg_nullcheck(argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  if (arg_parse (argc, argv, argtable)) {
    biomcmc_error ("If there are no other arguments, I can read everything from the control file.");
  }
  if (at_help->count) print_usage (argtable, argv[0], true);
  if (at_nsa1->count) ctrl->sa_n_samples  = at_nsa1->ival[0];
  if (at_lbda->count) {
    for (i=0; i < at_lbda->count; i++) ctrl->lambda_prior[i] = at_lbda->dval[i]; 
    for (; i < Ndist; i++) ctrl->lambda_prior[i] = at_lbda->dval[0]; 
  }
  arg_freetable (argtable, sizeof(argtable)/sizeof(argtable[0]));
  return ctrl;
}

void 
print_usage (void **argtable, char *progname)
{
  printf ("This program can use a control file with information about hyperprior parameters, gene file names,\n");
  printf ("  species names, etc. You can overwrite some or all of these parameters by using command line arguments,\n");
  printf ("  such that you can 'recycle' the control file for more than one analysis. One obvious example \n");
  printf ("  would be to run the MCMC and then analyse its output (using the same control file, with different\n");
  printf ("  command-line options). You can also define everything with command-line parameters.\n\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, argtable, "\n\n");
  arg_print_glossary(stdout, argtable,"  %-28s %s\n");
  exit (EXIT_FAILURE);
}

