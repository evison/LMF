
#include "metisbin.h"
#include "bmetis.h"

int main(int argc, char *argv[])
{
	idx_t options[METIS_NOPTIONS];
  	bigraph_t *bigraph;
  	idx_t *perm, *iperm;
  	params_t *params;
  	int status, i, j;

  	/* rdiags[i][0] and cdiags[i][0] saves the length of each array
  	 * excluding the first value */
  	idx_t **rdiags, **cdiags;
  	idx_t ndiags;

  	params = parse_cmdline(argc, argv);

  	gk_startcputimer(params->iotimer);
  	bigraph = ReadBiGraph(params);
  	gk_stopcputimer(params->iotimer);

  	if(bigraph == NULL){
  		printf("Input Error : nrows + ncols != nvtxs\n");
  		printf("\n***Metis returned with an error.\n");
  		return -1;
  	}

  	BDFPrintInfo(params, bigraph);

    METIS_SetDefaultOptions(options);
    /*User specific parameters*/
	options[METIS_OPTION_CTYPE]    = params->ctype;
	options[METIS_OPTION_IPTYPE]   = params->iptype;
	options[METIS_OPTION_RTYPE]    = params->rtype;
	options[METIS_OPTION_CCORDER]  = params->ccorder;
	options[METIS_OPTION_SEED]     = params->seed;
	options[METIS_OPTION_DBGLVL]   = params->dbglvl;
	options[METIS_OPTION_DENSITY] = params->density * DIVIDER;
	options[METIS_OPTION_NROWS] = params->nrows;
	options[METIS_OPTION_NCOLS] = params->ncols;
	options[METIS_OPTION_KAPPA] = params->kappa;
	options[METIS_OPTION_NDIAGS] = params->ndiags;

	/*Inner parameters*/
	options[METIS_OPTION_COMPRESS] = params->compress;
	options[METIS_OPTION_UFACTOR]  = params->ufactor;
	options[METIS_OPTION_PFACTOR]  = params->pfactor;
	options[METIS_OPTION_NCUTS] = params->ncuts;
	options[METIS_OPTION_NSEPS]    = params->nseps;
	options[METIS_OPTION_NITER]    = params->niter;
	options[METIS_OPTION_OBJTYPE] = params->objtype;

	perm  = imalloc(bigraph->super->nvtxs, "main: perm");
  	iperm = imalloc(bigraph->super->nvtxs, "main: iperm");

	gk_malloc_init();
	gk_startcputimer(params->parttimer);

	/* Initialize my global paramters */
	gk_clearcputimer(_parttimer);
	gk_clearcputimer(_nztimer);
	_totalcheck = 0;
	_firsthit = 0;
	_maxarea = -1;
	_maxnz = -1;
	_minarea = 700000000000;
	_minnz = 300000000;
	_avgarea = 0;
	_avgnz = 0;
	_maxdense = 0;
	_mindense = 0;

  	/* All the memory that is not allocated in this file should be allocated after
  	 * gk_malloc_init() and be freed before gk_GetCurMemoryUsed().
  	 * Memory that is allocated in this file should be free in the end of main()*/
  	status = METIS_NodeBDF(&bigraph->super->nvtxs, bigraph->super->xadj, bigraph->super->adjncy,
  			bigraph->super->vwgt, bigraph->nrows, bigraph->ncols,
  			options, bigraph->rlabel->label, bigraph->rlabel->ref, bigraph->clabel->label, bigraph->clabel->ref,
  			&rdiags, &cdiags, &ndiags, perm, iperm);

  	gk_stopcputimer(params->parttimer);

	if (gk_GetCurMemoryUsed() != 0)
    	printf("***It seems that Metis did not free all of its memory!\n");
	params->maxmemory = gk_GetMaxMemoryUsed();
	gk_malloc_cleanup(0);

	if (status != METIS_OK) {
		printf("\n***Metis returned with an error.\n");
	}
	else {
		if (! params->nooutput) {
	  		/* Write the permutation */
	  		gk_startcputimer(params->iotimer);
	  		WritePermutation(params->filename, iperm, bigraph->super->nvtxs);
	  		WriteDiags(params->filename, rdiags, cdiags, ndiags);
	  		gk_stopcputimer(params->iotimer);
		}
		BDFReportResults(params, bigraph);
	}

	/* free inner function memory */
	for (i = 0; i < ndiags; i++) {
		free((void*)rdiags[i]);
		free((void*)cdiags[i]);
	}
	free((void*)rdiags);
	free((void*)cdiags);

	/* free memroy allocated in this function */
	FreeBiGraph((ctrl_t*)NULL, &bigraph);
	gk_free((void **)&perm, &iperm, LTERM);
	gk_free((void **)&params->filename, &params->tpwgtsfile, &params->tpwgts,
	  &params->ubvec, &params, LTERM);

	return status;
}

/*************************************************************************/
/*! This function prints run parameters */
/*************************************************************************/
void BDFPrintInfo(params_t *params, bigraph_t *bigraph)
{
	printf("******************************************************************************\n");
	printf("%s", METISTITLE);
	printf(" (HEAD: %s, Built on: %s, %s)\n", SVNINFO, __DATE__, __TIME__);
	printf(" size of idx_t: %zubits, real_t: %zubits, idx_t *: %zubits\n",
		8*sizeof(idx_t), 8*sizeof(real_t), 8*sizeof(idx_t *));
	printf("\n");

	printf("Bipatite Graph Information --------------------------------------------------\n");
	printf(" Name: %s, #RowVertices: %"PRIDX", #ColumnVertices: %"PRIDX", #Edges: %"PRIDX"\n",
		params->filename, bigraph->nrows, bigraph->ncols, bigraph->super->nedges/2);

	printf("\n");
	printf("User Specific Options -------------------------------------------------------\n");
	printf(" ctype=%s, rtype=%s, iptype=%s, seed=%"PRIDX", dbglvl=%"PRIDX", ccorder=%s, compress=%s\n",
		ctypenames[params->ctype], rtypenames[params->rtype], iptypenames[params->iptype],
		params->seed, params->dbglvl, (params->ccorder  ? "YES" : "NO"), (params->compress ? "YES" : "NO"));
	printf(" density=%.4f, kappa=%d, nrows=%d, ncols=%d, area=%lld\n",
			params->density, params->kappa, params->nrows, params->ncols, params->nrows*params->ncols);

	printf("\n");
	printf("Inner Options ---------------------------------------------------------------\n");
	printf(" niter=%"PRIDX", nseps=%"PRIDX", ncuts=%"PRIDX", nparts=%"PRIDX"\n", params->niter,
			params->nseps, params->ncuts, params->nparts);
	printf(" ufactor=%.3f, pfactor=%.2f, nooutput=%s, objtype=%s\n\n", I2RUBFACTOR(params->ufactor),
			0.1*params->pfactor, (params->nooutput ? "YES" : "NO"), objtypenames[params->objtype] );
	printf("\n");
	printf("Node-based Nested Dissection ------------------------------------------------\n");
}


/*************************************************************************/
/*! This function does any post-ordering reporting */
/*************************************************************************/
void BDFReportResults(params_t *params, bigraph_t *bigraph)
{
	gk_startcputimer(params->reporttimer);
	gk_stopcputimer(params->reporttimer);

	printf("\nTiming Information ----------------------------------------------------------\n");
	printf("  I/O:          \t\t %7.3"PRREAL" sec\n", gk_getcputimer(params->iotimer));
	printf("  Ordering:     \t\t %7.3"PRREAL" sec   (METIS time)\n", gk_getcputimer(params->parttimer));
	printf("  Reporting:    \t\t %7.3"PRREAL" sec\n", gk_getcputimer(params->reporttimer));
	printf("  Partitioning: \t\t %7.3"PRREAL" sec\n", gk_getcputimer(_parttimer));
	printf("  NZStating:    \t\t %7.3"PRREAL" sec\n", gk_getcputimer(_nztimer));
	printf("\nMemory Information ----------------------------------------------------------\n");
	printf("  Max memory used:\t\t %7.3"PRREAL" MB\n", (real_t)(params->maxmemory/(1024.0*1024.0)));
	printf("\nHeuristic Information -------------------------------------------------------\n");
	printf("  TotalCheck:   \t\t %"PRIDX"\n", _totalcheck);
	printf("  FirstHit:     \t\t %"PRIDX"\n", _firsthit);
	printf("  FirstHitRate: \t\t %7.3"PRREAL"\n", (_totalcheck == 0 ? 1 : (real_t)1.0*_firsthit/_totalcheck));
	printf("  MaxArea:      \t\t %"PRIDX"\n", _maxarea);
	printf("  MaxNonZeros:  \t\t %"PRIDX"\n", _maxnz);
	printf("  MinArea:      \t\t %"PRIDX"\n", _minarea);
	printf("  MinNonZeors:  \t\t %"PRIDX"\n", _minnz);
	printf("  AvgArea:      \t\t %7.3"PRREAL"\n", _avgarea);
	printf("  AvgNz:        \t\t %7.3"PRREAL"\n", _avgnz);
	printf("  MaxDense:     \t\t %7.6"PRREAL"\n", _maxdense);
	printf("  MinDense:     \t\t %7.6"PRREAL"\n", _mindense);
	printf("******************************************************************************\n");

}
