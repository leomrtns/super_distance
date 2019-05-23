/* 
 * Copyright (C) 2019-today Leonardo de Oliveira Martins
 *
 * This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file src_super_distance.h
 *  \brief struct with argtable command line arguments 
 */

#ifndef _super_distance_src_control_file_h_
#define _super_distance_src_control_file_h_

#include <super_distance.h>   // these are the ../lib library headers 

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_lit  *fast;
  struct arg_dbl  *tol;
  struct arg_file *spname;
  struct arg_file *genfil;
  struct arg_file *outfil;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

#endif
