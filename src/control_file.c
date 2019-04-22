#include "control_file.h" 

//void print_usage (void **argtable, char *progname, struct arg_end *end);
void print_usage (arg_parameters params, char *progname);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .file = arg_filen("S","species" ,NULL, 1, argc-2, "file with species names, one name per line (nexus-style bracketed comments are allowed)"),
    .end  = arg_end(20)
  };
  void* argtable[] = {params.help, params.file, params.end};
  params.argtable = argtable;

  params.file->filename[0] = "out.txt";

  if (arg_nullcheck(argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  if (arg_parse (argc, argv, argtable)) print_usage (params, argv[0]);
  if (params.help->count) print_usage (params, argv[0]);
  
  return params;
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.end->count) {
    arg_print_errors(stdout, params.end, basename(progname));
    printf ("Error when reading arguments from command line\n\n");
  }
  printf ("This program can use a control file with information about hyperprior parameters, gene file names.\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-28s %s\n");
  arg_freetable (params.argtable, 3);
  exit (EXIT_FAILURE);
}

int
main (int argc, char **argv)
{
  arg_parameters params = get_parameters_from_argv (argc, argv);
  printf ("result: %s\n", params.file->filename[0]);
  arg_freetable (params.argtable, 3);
  return EXIT_SUCCESS;
}
