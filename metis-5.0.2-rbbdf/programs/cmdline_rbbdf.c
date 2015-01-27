/*!
\file cmdline_ndmetis.c
\brief Command-line argument parsing for ndmetis

\date 12/24/2008
\author George
\version\verbatim $Id: cmdline_ndmetis.c 10481 2011-07-05 18:01:23Z karypis $\endverbatim
*/

#include "metisbin.h"


/*-------------------------------------------------------------------
 * Command-line options
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
  {"ctype",          1,      0,      METIS_OPTION_CTYPE},
  {"iptype",         1,      0,      METIS_OPTION_IPTYPE},
  {"rtype",          1,      0,      METIS_OPTION_RTYPE},
  {"ufactor",        1,      0,      METIS_OPTION_UFACTOR},
  {"pfactor",        1,      0,      METIS_OPTION_PFACTOR},
  {"nocompress",     0,      0,      METIS_OPTION_COMPRESS},
  {"ccorder",        0,      0,      METIS_OPTION_CCORDER},
  {"nooutput",       0,      0,      METIS_OPTION_NOOUTPUT},
  {"niter",          1,      0,      METIS_OPTION_NITER},
  {"nseps",          1,      0,      METIS_OPTION_NSEPS},
  {"seed",           1,      0,      METIS_OPTION_SEED},
  {"dbglvl",         1,      0,      METIS_OPTION_DBGLVL},
  {"help",           0,      0,      METIS_OPTION_HELP},
  {"ncuts",          1,      0,      METIS_OPTION_NCUTS},
  {"objtype",        1,      0,      METIS_OPTION_OBJTYPE},
	/*evison*/
  {"density",        1,      0,      METIS_OPTION_DENSITY},	//1 means that this parameter needs input
  {"kappa",          1,      0,      METIS_OPTION_KAPPA},
  {"nrows",          1,      0,      METIS_OPTION_NROWS},
  {"ncols",          1,      0,      METIS_OPTION_NCOLS},
  {"ndiags",         1,      0,      METIS_OPTION_NDIAGS},
  {0,                0,      0,      0}
};


static gk_StringMap_t ctype_options[] = {
 {"rm",                 METIS_CTYPE_RM},
 {"shem",               METIS_CTYPE_SHEM},
 {NULL,                 0}
};

static gk_StringMap_t iptype_options[] = {
 {"edge",               METIS_IPTYPE_EDGE},
 {"node",               METIS_IPTYPE_NODE},
 {"grow",               METIS_IPTYPE_GROW},
 {NULL,                 0}
};

static gk_StringMap_t rtype_options[] = {
 {"2sided",             METIS_RTYPE_SEP2SIDED},
 {"1sided",             METIS_RTYPE_SEP1SIDED},
 {"fm",                 METIS_RTYPE_FM},
 {NULL,                 0}
};

static gk_StringMap_t objtype_options[] = {
 {"cut",                METIS_OBJTYPE_CUT},
 {"vol",                METIS_OBJTYPE_VOL},
 {"node",               METIS_OBJTYPE_NODE},
 {NULL,                 0}
};


/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: rbbdf [options] <filename>",
" ",
" Required parameters",
"    filename    Stores the graph to be partitioned.",
" ",
" Optional parameters",
"  -density",
"     The density requirment.",
" ",
"  -nrows",
"     Number of rows vertices in the bipartite graph",
" ",
"  -ncols",
"     Number of column vertices in the bipartite graph",
" ",
"  -kappa",
"     Number of vertex separators to compute in each round,",
"     the best one will be taken as the final separator",
" ",
"  -ccorder",
"     Extract connected components before separation",
" ",
"  -iptype=string [applies only when -ptype=rb]",
"     Specifies the scheme to be used to compute the initial bisection",
"     of the graph.",
"     The possible values are:",
"        edge     - Separator from an edge cut",
"        node     - Separator from a greedy node-based strategy [default]",
" ",
"  -ctype=string",
"     Specifies the scheme to be used to match the vertices of the graph",
"     during the coarsening.",
"     The possible values are:",
"        rm       - Random matching",
"        shem     - Sorted heavy-edge matching [default]",
" ",
"  -rtype=string",
"     Specifies the scheme to be used for refinement.",
"     The possible values are:",
"        1sided   - 1-sided node-based refinement [default]",
"        2sided   - 2-sided node-based refinement",
" ",
"  -seed=int",
"     Selects the seed of the random number generator.  ",
" ",
"  -dbglvl=int      ",
"     Selects the dbglvl.  ",
" ",
"  -help",
"     Prints this message.",
""
/*
"  -nocompress",
"     Do not compress a graph and do separation directly",
" ",
*/
};

static char shorthelpstr[][100] = {
" ",
"   Usage: rbbdf [options] <filename>",
"          use 'rbbdf -help' for a summary of the options.",
""
};

/*************************************************************************/
/*! This is the entry point of the command-line argument parser */
/*************************************************************************/
params_t *parse_cmdline(int argc, char *argv[])
{
  int i, j, k;
  int c, option_index;
  params_t *params;

  params = (params_t *)gk_malloc(sizeof(params_t), "parse_cmdline");
  memset((void *)params, 0, sizeof(params_t));

  /*fixed parameters*/
  params->ufactor       = OMETIS_DEFAULT_UFACTOR;
  params->pfactor       = 0;
  params->nooutput      = 0;
  params->wgtflag       = 1;
  params->nseps         = 1;
  params->niter         = 10;	/* evison */ /* XXX magic number */
  params->nparts        = 3;	/* evison */ /* TODO Note, the parameter is really set in METIS_NodeBDF() */
  /*params->ncuts = 10;*/		/* evison */ /* ncuts is not used */
  params->objtype       = METIS_OBJTYPE_NODE;
  params->filename      = NULL;
  params->compress      = 0;

  /*user specific parameters*/
  params->ctype         = METIS_CTYPE_SHEM;
  params->iptype        = METIS_IPTYPE_NODE;
  params->rtype         = METIS_RTYPE_SEP1SIDED;
  params->ccorder       = 0;
  params->seed          = -1;
  params->dbglvl        = 0;
  /* evison */
  params->density = -1;
  params->nrows = -1;
  params->ncols = -1;
  params->kappa = 1;
  params->ndiags = -1;

  gk_clearcputimer(params->iotimer);
  gk_clearcputimer(params->parttimer);
  gk_clearcputimer(params->reporttimer);

  /* Parse the command line arguments  */
  while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
      case METIS_OPTION_CTYPE:
        if (gk_optarg)
          if ((params->ctype = gk_GetStringID(ctype_options, gk_optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, gk_optarg);
        break;
      case METIS_OPTION_IPTYPE:
        if (gk_optarg)
          if ((params->iptype = gk_GetStringID(iptype_options, gk_optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, gk_optarg);
        break;

      case METIS_OPTION_RTYPE:
        if (gk_optarg)
          if ((params->rtype = gk_GetStringID(rtype_options, gk_optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, gk_optarg);
        break;
/*
      case METIS_OPTION_COMPRESS:
        params->compress = 0;
        break;
*/
      case METIS_OPTION_CCORDER:
        params->ccorder = 1;
        break;

      case METIS_OPTION_SEED:
        if (gk_optarg) params->seed = (idx_t)atoi(gk_optarg);
        break;

      case METIS_OPTION_DBGLVL:
        if (gk_optarg) params->dbglvl = (idx_t)atoi(gk_optarg);
        break;

        /*evison*/
      case METIS_OPTION_DENSITY:
    	  if (gk_optarg) params->density = (real_t)atof(gk_optarg);
    	  break;

      case METIS_OPTION_KAPPA:
    	  if (gk_optarg) params->kappa = (idx_t)atoi(gk_optarg);
    	  break;

      case METIS_OPTION_NROWS:
    	  if (gk_optarg) params->nrows = (idx_t)atoi(gk_optarg);
    	  break;

      case METIS_OPTION_NCOLS:
    	  if (gk_optarg) params->ncols = (idx_t)atoi(gk_optarg);
    	  break;

      case METIS_OPTION_NDIAGS:
    	  if (gk_optarg) params->ndiags = (idx_t)atoi(gk_optarg);
    	  break;

      case METIS_OPTION_HELP:
        for (i=0; strlen(helpstr[i]) > 0; i++)
          printf("%s\n", helpstr[i]);
        exit(0);
        break;
      case '?':
      default:
        errexit("Illegal command-line option(s)\n"
                "Use %s -help for a summary of the options.\n", argv[0]);
    }
  }

  if (argc-gk_optind != 1) {
    printf("Missing parameters.");
    for (i=0; strlen(shorthelpstr[i]) > 0; i++)
      printf("%s\n", shorthelpstr[i]);
    exit(0);
  }

  params->filename = gk_strdup(argv[gk_optind++]);

  return params;
}


