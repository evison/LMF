#include "metislib.h"
#include "bmetis.h"
#include "sys/malloc.h"
/*************************************************************************/
/*! This function is the entry point for the multilevel nested dissection
    ordering code. At each bisection, a node-separator is computed using
    a node-based refinement approach.

    \param nvtxs is the number of vertices in the graph.
    \param xadj is of length nvtxs+1 marking the start of the adjancy
           list of each vertex in adjncy.
    \param adjncy stores the adjacency lists of the vertices. The adjnacy
           list of a vertex should not contain the vertex itself.
    \param vwgt is an array of size nvtxs storing the weight of each
           vertex. If vwgt is NULL, then the vertices are considered
           to have unit weight.
    \param numflag is either 0 or 1 indicating that the numbering of
           the vertices starts from 0 or 1, respectively.
    \param options is an array of size METIS_NOPTIONS used to pass
           various options impacting the of the algorithm. A NULL
           value indicates use of default options.
    \param perm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A'[i] = A[perm[i]].
    \param iperm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A[i] = A'[iperm[i]].
*/
/*************************************************************************/
int METIS_NodeBDF(idx_t *nvtxs, idx_t* xadj, idx_t* adjncy, idx_t *vwgt, idx_t nrows, idx_t ncols,
		idx_t *options, idx_t *rlabel_label, idx_t rlabel_ref, idx_t *clabel_label, idx_t clabel_ref,
		idx_t ***r_rdiags, idx_t ***r_cdiags, idx_t *r_ndiags, idx_t *perm, idx_t *iperm){

	int sigrval = 0, renumber=0;
	ctrl_t *ctrl;
	idx_t *cptr, *cind;
	bigraph_t *obigraph = NULL;
	graph_t *ograph = NULL;
	idx_t nnvtxs;
	label_t *rlabel, *clabel;
	int i, j;

	/* set up malloc cleaning code and signal catchers */
	if (!gk_malloc_init())	return METIS_ERROR_MEMORY;

	gk_sigtrap();

	if ((sigrval = gk_sigcatch()) != 0)	goto SIGTHROW;

	/*1-ncon 3-nparts*/	/*TODO can nparts be assigned other numbers?*/
	ctrl = SetupCtrl(METIS_OP_BMETIS, options, 1, 3, NULL, NULL);
	if (!ctrl) {
		gk_siguntrap();
		gk_malloc_cleanup(0);
		return METIS_ERROR_INPUT;
	}

	/* if required, change the numbering to 0 */
	if (ctrl->numflag == 1) {
		Change2CNumbering(*nvtxs, xadj, adjncy);
		renumber = 1;
	}

	IFSET(ctrl->dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
  	IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->TotalTmr));

  	/* prune the dense columns */
	if (ctrl->pfactor > 0.0) {	/* TODO */
		printf("**Waring : Dense columns pruning should not be performed\n");
		ctrl->pfactor = 0.0;
	}

	rlabel = (label_t*)gk_malloc(sizeof(label_t), "METIS_NodeBDF : rlabel");
	clabel = (label_t*)gk_malloc(sizeof(label_t), "METIS_NodeBDF : clabel");
	rlabel->label = rlabel_label;	rlabel->ref = rlabel_ref;
	clabel->label = clabel_label;	clabel->ref = clabel_ref;
	ograph = SetupGraph(ctrl, *nvtxs, 1, xadj, adjncy, vwgt, NULL, NULL);

	obigraph = SetupBiGraphFromGraph(ograph, nrows, ncols, ctrl->nrows, ctrl->ncols, rlabel, clabel);
	ctrl->obigraph = obigraph;

	ASSERT(CheckBiGraph(obigraph, ctrl->numflag, 1));	/* TODO debug */

	/* allocate workspace memory */
  	AllocateWorkSpace(ctrl, obigraph->super);

  	if (ctrl->ccorder)
		MlevelNestedBDFCC(ctrl, obigraph, iperm, 1, r_rdiags, r_cdiags, r_ndiags);
	else
		MlevelNestedBDF(ctrl, obigraph, iperm, 1, r_rdiags, r_cdiags, r_ndiags);

	for (i = 0; i < *nvtxs; i++)
		perm[iperm[i]] = i;

	ASSERT(CheckPermIPerm(perm, iperm, *nvtxs));	/* TODO debug */

	IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->TotalTmr));
	IFSET(ctrl->dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

	/* clean up */
	FreeBiGraph(ctrl, &ctrl->obigraph);
	FreeCtrl(&ctrl);

SIGTHROW:
	/* if required, change the numbering back to 1 */
	if (renumber) Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);

	gk_siguntrap();
	gk_malloc_cleanup(0);

	return metis_rcode(sigrval);
}


/***********************************************************************
 * This functhion performs real permutation and writes resutls to
 * data structures given by parameters.
 * 	\param ctrl is the ctrl parameters for reordering
 * 	\param head is head of block list to reorder
 * 	\param order is for storing and returning current order of vertices
 * 	\param lastvtx is the last vertex index in this round
 * 	\param rdiags is for returning the row indices of block diagonals
 * 	\param cdiags is for returning the row indices of block diagonals
 * 	\param ndiags is for returning number of block diagonal extracted
 *
 *	TODO In order to save memory, bigraph can be deleted in this
 *	function if a patition has been performed and the graph is split
 *	into several block diagonals, because we only need to keep track
 *	of current diagonals. But if so, we must make sure the memory is
 *	not refreed in main() function.
 *
 *	bigraph, rdiags, cdiags, ndiags have been initialized well.
 ***********************************************************************/
void MlevelNestedBDF(ctrl_t *ctrl, bigraph_t *head, idx_t *order, idx_t ndiags,
		idx_t ***r_rdiags, idx_t ***r_cdiags, idx_t *r_ndiags){

	real_t avgdensity, replacedensity, mindensity, tempdensity;
	bigraph_t *lbigraph, *rbigraph, *temp, *p, *q;
	bigraph_t **sort;
	real_t *denses;
	idx_t *areas;
	graph_t *cgraph, *tosplit;
	idx_t *cptr, *cind, *bndind, *label;
	idx_t i, j, k, pos, nnvtxs, nbnd;
	graph_t *lgraph, *rgraph, *swap;
	idx_t partitioned = 0;

	avgdensity = AverageDensity(head);
	printf("ndiags = %d, avgdensity = %.6f, requirdensity = %.6f\n", ndiags, avgdensity, ctrl->density);

	if( (ctrl->ndiags != -1 && ndiags >= ctrl->ndiags)
			|| (ctrl->ndiags == -1 && avgdensity >= ctrl->density) ){
		/* TODO debug */
		sort = (bigraph_t**)gk_malloc(ndiags * sizeof(bigraph_t*), "MlevelNestedBDF : sorted diagonalas");

		areas = (idx_t*)gk_malloc(ndiags * sizeof(idx_t), "MlevelNestedBDF : area list");

		SortBlockDiagsBySingleArea(head, ndiags, sort, areas);

		OrderEachGraph(head, order);
		ConstructResult(head, ndiags, r_rdiags, r_cdiags, r_ndiags);

		printf("***RETURN @1: Density requiment reached\n");
		return;
	}

	/* Sort the block diagonals according to their densities, and save the sorted densities
	 * in 'densities' and their corresponding pointers in 'sort'
	 * TODO which block diagonal to choose? */

	sort = (bigraph_t**)gk_malloc(ndiags * sizeof(bigraph_t*), "MlevelNestedBDF : sorted diagonalas");

	areas = (idx_t*)gk_malloc(ndiags * sizeof(idx_t), "MlevelNestedBDF : area list");

	SortBlockDiagsBySingleArea(head, ndiags, sort, areas);

	_totalcheck++;	/* TODO exp heuristic */
	/* try to split the block diagonals from the one with least average density on,
	 * until the average density of the resulting blocks is improved*/
	for (i = 0; i < ndiags; i++) {
		if (!sort[i]->partible) continue;
		/* The graph could not be compressed on the whole bipartie graph in advance,
		 * because if an row vetex and a column vertex are compressed into one vertex,
		 * then this would not be a biparatite graph any more !
		 * Compress could be done here before each real partitioning, and the partition
		 * result should be mapped back into orignal indices once finished. */
		if (ctrl->compress) {
			cptr = imalloc(sort[i]->super->nvtxs+1, "BMETIS: cptr");
			cind = imalloc(sort[i]->super->nvtxs, "BMETIS: cind");

			cgraph = CompressGraph(ctrl, sort[i]->super->nvtxs, sort[i]->super->xadj,
					sort[i]->super->adjncy, sort[i]->super->vwgt, cptr, cind);
			if (cgraph == NULL) {
				/* if there was no compression, cleanup the compressed flag */
				gk_free((void **)&cptr, &cind, LTERM);
				ctrl->compressed = 0;
				tosplit = sort[i]->super;
			}
			else {
				nnvtxs = sort[i]->super->nvtxs;
				ctrl->compressed = 1;
				tosplit = cgraph;
				/*
				ctrl->cfactor = 1.0*(*graph->nvtxs)/nnvtxs;
				if (ctrl->cfactor > 1.5 && ctrl->nseps == 1)	ctrl->nseps = 2;
				*//* TODO */
			}
		}
		else {
			ctrl->compressed = 0;
			tosplit = sort[i]->super;
		}

		ASSERT(CheckGraph(tosplit, ctrl->numflag, 1));	/* TODO debug */

		gk_startcputimer(_parttimer);	/* TODO debug timer */
		MlevelNodeBisectionMultipleBDF(ctrl, tosplit);
		gk_stopcputimer(_parttimer);

		IFSET(ctrl->dbglvl, METIS_DBG_SEPINFO,
      		printf("Nvtxs: %6"PRIDX", [%6"PRIDX" %6"PRIDX" %6"PRIDX"]\n",
        		sort[i]->super->nvtxs, sort[i]->super->pwgts[0],
        		sort[i]->super->pwgts[1], sort[i]->super->pwgts[2]));

		/* extract resulting blocks diagonal list, note that compression may have been done */
		if (ctrl->compressed) {
			SplitGraphOrderUncompressBDF(ctrl, sort[i]->super, cgraph, cptr, cind, &lgraph, &rgraph);
			FreeGraph(&cgraph);
			gk_free((void **)&cptr, &cind, LTERM);
		}
		else {
			gk_startcputimer(_parttimer);	/* TODO exp timer */
			SplitGraphOrderBDF(ctrl, sort[i]->super, &lgraph, &rgraph);
			/*swap = lgraph; lgraph = rgraph; rgraph = swap;*/
			gk_stopcputimer(_parttimer);
		}

		if (lgraph->nvtxs == 0 || rgraph->nvtxs == 0){
			sort[i]->partible = 0;
			FreeGraph(&lgraph);
			FreeGraph(&rgraph);
			continue;
		}

		/* construct new left and right bigraphs */
		ExtractBiGraph(ctrl, sort[i], lgraph, rgraph, &lbigraph, &rbigraph);

		/* check whether average density is improved */
		replacedensity = AverageReplaceDensity (head, sort[i], lbigraph, rbigraph);

		if (replacedensity > avgdensity) {	/* if so */
			if (i == 0)	_firsthit++;	/* TODO exp heuristic */
			partitioned = 1;
			/* insert new bigraphs into block diagonal list */
			p = head;	while (p && p != sort[i])	p = p->next;
			if (p == NULL) { printf("***ERROR: p == NULL"); exit(-1); }	/* TODO debug */
			if (p == head) { /* p is head */
				head = lbigraph;
				lbigraph->next = rbigraph;
				rbigraph->next = p->next;
			}
			else {
				q = head;	while (q->next != p)	q = q->next;
				q->next = lbigraph;
				lbigraph->next = rbigraph;
				rbigraph->next = p->next;
			}
			for (j = 0; j < sort[i]->super->nbnd; j++) {	/* manage order */
				order[sort[i]->super->label[sort[i]->super->bndind[j]]] = --sort[i]->lastvtx;
			}
			rbigraph->lastvtx = sort[i]->lastvtx;
			lbigraph->lastvtx = rbigraph->lastvtx - rbigraph->super->nvtxs;

			/* release sort[i], but ctrl->obigraph should not be released */
			FreeBiGraph(ctrl, &(sort[i]));
			break;
		}
		else {	/* else */
			/* release memory of lbigrahp and rbigraph then try the next block diagonal */
			FreeBiGraph(ctrl, &lbigraph);
			FreeBiGraph(ctrl, &rbigraph);
		}
	}

	gk_free((void**)&sort, &areas, LTERM);

	if (partitioned) { /* recursive call */
		MlevelNestedBDF(ctrl, head, order, ndiags+1, r_rdiags, r_cdiags, r_ndiags);
	}
	else {	/* non of the diagonal block improves average density */
		/* manage order of each block diagonal graph */
		OrderEachGraph(head, order);
		ConstructResult(head, ndiags, r_rdiags, r_cdiags, r_ndiags);

		printf("***RETURN @2: Can not improve density any more.\n");
		return;
	}
}

/***********************************************************************
 * This functhion is similar to its non 'CC' version execpt that it
 * extracts connected components of a graph before reording.
 * 	\param ctrl is the ctrl parameters for reordering
 * 	\param head is the head of block list to reorder
 * 	\param order is for storing and returning current order of vertices
 * 	\param lastvtx is the last vertex index in this round
 * 	\param rdiags is for returning the row indices of block diagonals
 * 	\param cdiags is for returning the row indices of block diagonals
 * 	\param ndiags is for returning number of block diagonal extracted
 *
 ***********************************************************************/
void MlevelNestedBDFCC(ctrl_t *ctrl, bigraph_t *head, idx_t *order, idx_t ndiags,
		idx_t *** r_rdiags, idx_t ***r_cdiags, idx_t *r_ndiags){

}

/**
 * This function orders the permutation for each graph in order.
 * Note, all the borders nodes have been ordered well when a
 * partition was taken, so we don't need to permute border nodes
 * in this function.
 */
void OrderEachGraph(bigraph_t *head, idx_t *order) {
	idx_t i;
	bigraph_t *p = head;

	while (p) {
		for (i = 0; i < p->super->nvtxs; i++) {
			order[p->super->label[i]] = --p->lastvtx;
		}
		p = p->next;
	}
}

void ConstructResult(bigraph_t *head, idx_t ndiags, idx_t ***r_rdiags, idx_t ***r_cdiags, idx_t *r_ndiags) {
	idx_t **rdiags, **cdiags;
	bigraph_t *p, *q;
	idx_t i, j, k, nrows, ncols, snrows, sncols;

	idx_t sarea = 0, snz = 0, area, nz, cnt = 0;
	real_t dense;

	/* used malloc instead of gk_malloc, becasue gk_malloced memory will be released
	 * in gkmalloc_clean_up in METIS_NodeBDF, but we stil nedd the memory in rbbdf.c */
	rdiags = (idx_t**)malloc(ndiags * sizeof(idx_t*));
	cdiags = (idx_t**)malloc(ndiags * sizeof(idx_t*));

	i = 0;
	p = head;
	while (p) {
		nrows = 0;	q = p;
		while (q) {
			nrows += q->nrows;
			q = q->down;
		}

		ncols = 0;	q = p;
		while (q) {
			ncols += q->ncols;
			q = q->right;
		}

		rdiags[i] = (idx_t*)malloc((nrows + 1) * sizeof(idx_t));
		cdiags[i] = (idx_t*)malloc((ncols + 1) * sizeof(idx_t));
		rdiags[i][0] = nrows;
		cdiags[i][0] = ncols;

		j = 1;	q = p;
		while (q) {
			for (k = 0; k < q->nrows; k++)	rdiags[i][j++] = q->rlabel->label[k];
			q = q->down;
		}

		j = 1;	q = p;
		while (q) {
			for (k = 0; k < q->ncols; k++)	cdiags[i][j++] = q->clabel->label[k];
			q = q->right;
		}

		p = p->next;
		i++;
	}

	*r_ndiags = ndiags;
	*r_rdiags = rdiags;
	*r_cdiags = cdiags;

	/* Some statistics */
	p = head;
	i = cnt = snrows = sncols = 0;
	while (p) {
		cnt++;
		StatNzAndArea(p, &nz, &area, 0);
		sarea += area;
		snz += nz;
		dense = 1.0 * nz / area;

		if (area > _maxarea)	_maxarea = area;
		if (nz > _maxnz)	_maxnz = nz;
		if (area < _minarea)	_minarea = area;
		if (nz < _minnz)	_minnz = nz;
		if (dense > _maxdense)	_maxdense = dense;
		if (dense < _mindense)	_mindense = dense;

		StatNrowsAndNcols(p, &nrows, &ncols);
		snrows += nrows;
		sncols += ncols;
		p = p->next;
	}

	printf("Average nrows = %.2f, Average ncols = %.2f\n", 1.0*snrows/cnt, 1.0*sncols/cnt);

	_avgarea = 1.0 * sarea / cnt;
	_avgnz = 1.0 * snz / cnt;
}

/**
 * This function calculates the total area and non-zeros of a graph.
 * It DOES take borders into account.
 */
void StatNzAndArea(bigraph_t *bigraph, idx_t *r_snz, idx_t *r_sarea, idx_t islist) {
	bigraph_t *pblock, *prow, *pcolumn;
	idx_t snz = 0, sarea = 0;
	idx_t nrblocks = 0, cnt = 0;

	pblock = bigraph;
	while(pblock){
		cnt = 0;
		nrblocks = 0;
		prow = pblock;
		while(prow){
			nrblocks++;
			pcolumn = prow;
			while(pcolumn){
				cnt++;
				snz += pcolumn->nz;
				sarea += pcolumn->area;
				pcolumn = pcolumn->right;
			}
			prow = prow->down;
		}
		if (! islist)	break;
		pblock = pblock->next;

		if (cnt != nrblocks * nrblocks) {	/* TODO debug */
			printf("***ERROR: # blocks is not a square number !\n"); exit(-1);
		}
	}

	*r_snz = snz;
	*r_sarea = sarea;
}

void StatNrowsAndNcols(bigraph_t *bigraph, idx_t *r_nrows, idx_t *r_ncols) {
	idx_t nrows = 0, ncols = 0;
	bigraph_t *p;

	p = bigraph;
	while (p) {
		nrows += p->nrows;
		p = p->down;
	}

	p = bigraph;
	while (p) {
		ncols += p->ncols;
		p = p->right;
	}

	*r_nrows = nrows;
	*r_ncols = ncols;
}

/**
 * This functions computes the average density of block diagonals.
 * Note that it DOES take borders into account.
 */
real_t AverageDensity(bigraph_t *head) {
	idx_t snz = 0, sarea = 0;
	StatNzAndArea(head, &snz, &sarea, 1);
	return 1.0 * snz / sarea;
}

/**
 * This function calculates the new average denstiy if we replace bigraph old with
 * new bigraphs new1 and new2.
 */
real_t AverageReplaceDensity(bigraph_t *head, bigraph_t *old, bigraph_t *new1, bigraph_t *new2) {
	idx_t snz = 0, sarea = 0, nz = 0, area = 0;

	StatNzAndArea(head, &nz, &area, 1);
	snz += nz;	sarea += area;

	StatNzAndArea(old, &nz, &area, 0);
	snz -= nz;	sarea -= area;

	StatNzAndArea(new1, &nz, &area, 0);
	snz += nz;	sarea += area;

	StatNzAndArea(new2, &nz, &area, 0);
	snz += nz;	sarea += area;

	return 1.0 * snz / sarea;
}

/*
 * This function computes density of a bipartitie graph.
 * Note that it DOES take borders into account.
 */
real_t Density(bigraph_t *bigraph)
{
	idx_t snz = 0, sarea = 0;
	StatNzAndArea(bigraph, &snz, &sarea, 0);

	return 1.0 * snz / sarea;
//	return 1.0 * bigraph->nz / bigraph->area;
}

/*
 * This function computes area of a bipartitie graph.
 * Note that it DOES take borders into account.
 */
idx_t Area(bigraph_t *bigraph)
{
	bigraph_t *p, *q;
	idx_t sarea;

	sarea = 0;	p = bigraph;
	while (p) {
		q = p;
		while (q) {
			sarea += q->area;
			q = q->right;
		}
		p = p->down;
	}

	return sarea;
}

/*
 * This functions returns a sorted list of bigraph_t* pointers
 * according to the density of each block diagonal.
 */
void SortBlockDiagsByDense(bigraph_t *head, idx_t ndiags, bigraph_t** sort, real_t* denses){
	bigraph_t *pblock = head;
	idx_t i = 0;

	while(pblock){
		sort[i] = pblock;
		denses[i] = Density(pblock);
		i++;
		pblock = pblock->next;
	}

	ASSERT(i == ndiags);

	QSortBlockDiagsByDense(sort, denses, 0, ndiags-1);
}

/*
 * This function quick sorts 'densities', and at the same time arrange the pointer
 * list 'sort' according to 'densiteis', from the smallest to the largest
 */
void QSortBlockDiagsByDense(bigraph_t** sort, real_t* denses, idx_t l, idx_t r) {
    idx_t i, j;
    real_t x;
    bigraph_t *xx;

    if (l < r) {
        i = l;
        j = r;
        x = denses[i];
        xx = sort[i];
        while (i < j) {
            while(i < j && denses[j] > x) j--;
            if(i < j) {
            	denses[i] = denses[j];
				sort[i] = sort[j];
            	i++;
            }
            while(i < j && denses[i] < x) i++;
            if(i < j){
            	denses[j] = denses[i];
            	sort[j] = sort[i];
            	j--;
            }
        }
        denses[i] = x;
        sort[i] = xx;
        QSortBlockDiagsByDense(sort, denses, l, i-1);
        QSortBlockDiagsByDense(sort, denses, i+1, r);
    }
}

/*
 * This functions returns a sorted list of bigraph_t* pointers
 * according to the area of each block diagonal, namely, it DOES take
 * bordres into account.
 */
void SortBlockDiagsByArea(bigraph_t *head, idx_t ndiags, bigraph_t **sort, idx_t *areas){
	bigraph_t *pblock = head;
	idx_t i = 0;

	while(pblock){
		sort[i] = pblock;
		areas[i] = Area(pblock);
		i++;
		pblock = pblock->next;
	}

	ASSERT(i == ndiags);

	QSortBlockDiagsByArea(sort, areas, 0, ndiags-1);
}

/**
 * This function resturns a sorted list of bigraph_t* pointers
 * according to the single area of each block diagnal, namely, it DOES NOT take
 * borders into account.
 */
void SortBlockDiagsBySingleArea(bigraph_t *head, idx_t ndiags, bigraph_t **sort, idx_t *areas){
	bigraph_t *pblock = head;
	idx_t i = 0;

	while(pblock){
		sort[i] = pblock;
		areas[i] = pblock->area;
		i++;
		pblock = pblock->next;
	}

	ASSERT(i == ndiags);

	QSortBlockDiagsByArea(sort, areas, 0, ndiags-1);
}

/*
 * This function quick sorts 'areas', and at the same time arrange the pointer
 * list 'sort' according to 'areas', from the largest to the smallest
 */
void QSortBlockDiagsByArea(bigraph_t** sort, idx_t *areas, idx_t l, idx_t r) {
    idx_t i, j, x;
    bigraph_t *xx;

    if (l < r) {
        i = l;
        j = r;
        x = areas[i];
        xx = sort[i];
        while (i < j) {
            while(i < j && areas[j] < x) j--;
            if(i < j) {
            	areas[i] = areas[j];
				sort[i] = sort[j];
            	i++;
            }
            while(i < j && areas[i] > x) i++;
            if(i < j){
            	areas[j] = areas[i];
            	sort[j] = sort[i];
            	j--;
            }
        }
        areas[i] = x;
        sort[i] = xx;
        QSortBlockDiagsByArea(sort, areas, l, i-1);
        QSortBlockDiagsByArea(sort, areas, i+1, r);
    }
}

/*************************************************************************/
/*! This function performs multilevel node bisection (i.e., tri-section).
    It performs multiple bisections and selects the best. */
/*************************************************************************/
void MlevelNodeBisectionMultipleBDF(ctrl_t *ctrl, graph_t *graph)
{
	idx_t i, mincut;
	idx_t *bestwhere;

	/* if the graph is small, just find a single vertex separator */
	if (ctrl->nseps == 1 || graph->nvtxs < (ctrl->compress ? 1000 : 2000)) {
		MlevelNodeBisectionBDFL2(ctrl, graph, LARGENIPARTS);
		return;
	}

	WCOREPUSH;

	bestwhere = iwspacemalloc(ctrl, graph->nvtxs);

	mincut = graph->tvwgt[0];
	for (i=0; i<ctrl->nseps; i++) {
		MlevelNodeBisectionBDFL2(ctrl, graph, LARGENIPARTS);

		if (i == 0 || graph->mincut < mincut) {
			mincut = graph->mincut;
			if (i < ctrl->nseps-1)
				icopy(graph->nvtxs, graph->where, bestwhere);
		}

		if (mincut == 0)	break;

		if (i < ctrl->nseps-1)
			FreeRData(graph);
	}

	if (mincut != graph->mincut) {
		icopy(graph->nvtxs, bestwhere, graph->where);
		Compute2WayNodePartitionParams(ctrl, graph);
	}

	WCOREPOP;
}

/*************************************************************************/
/*! This function performs multilevel node bisection (i.e., tri-section).
    It performs multiple bisections and selects the best. */
/*************************************************************************/
void MlevelNodeBisectionBDFL2(ctrl_t *ctrl, graph_t *graph, idx_t niparts)
{
	idx_t i, mincut, nruns = 2;	/* XXX magic number! */
	graph_t *cgraph;
	idx_t *bestwhere;

	/* if the graph is small, just find a single vertex separator */
	if (graph->nvtxs < 5000) {
		MlevelNodeBisectionBDFL1(ctrl, graph, niparts);
		return;
	}

	WCOREPUSH;

	ctrl->CoarsenTo = gk_max(100, graph->nvtxs/30);
	cgraph = CoarsenGraphNlevels(ctrl, graph, 4);	/* XXX magic number! */

	bestwhere = iwspacemalloc(ctrl, cgraph->nvtxs);

	mincut = graph->tvwgt[0];
	for (i=0; i<nruns; i++) {
		MlevelNodeBisectionBDFL1(ctrl, cgraph, 0.7*niparts);
		if (i == 0 || cgraph->mincut < mincut) {
			mincut = cgraph->mincut;
			if (i < nruns-1)	icopy(cgraph->nvtxs, cgraph->where, bestwhere);
		}
		if (mincut == 0)	break;
		if (i < nruns-1)	FreeRData(cgraph);
	}

	if (mincut != cgraph->mincut)
	icopy(cgraph->nvtxs, bestwhere, cgraph->where);

	WCOREPOP;

	Refine2WayNode(ctrl, graph, cgraph);
}

/*************************************************************************/
/*! The top-level routine of the actual multilevel node bisection */
/*************************************************************************/
void MlevelNodeBisectionBDFL1(ctrl_t *ctrl, graph_t *graph, idx_t niparts)
{
	graph_t *cgraph;

	ctrl->CoarsenTo = graph->nvtxs/8;
	if (ctrl->CoarsenTo > 100)
	ctrl->CoarsenTo = 100;
	else if (ctrl->CoarsenTo < 40)
	ctrl->CoarsenTo = 40;

	cgraph = CoarsenGraph(ctrl, graph);

	niparts = gk_max(1, (cgraph->nvtxs <= ctrl->CoarsenTo ? niparts/2: niparts));

	InitSeparator(ctrl, cgraph, niparts);

	Refine2WayNode(ctrl, graph, cgraph);
}

/*************************************************************************/
/*! This function takes a graph and a tri-section (left, right, separator)
    and splits it into two graphs.

    This function relies on the fact that adjwgt is all equal to 1.
*/
/*************************************************************************/
void SplitGraphOrderBDF(ctrl_t *ctrl, graph_t *graph, graph_t **r_lgraph, graph_t **r_rgraph)
{
	idx_t i, ii, j, k, l, istart, iend, mypart, nvtxs, snvtxs[3], snedges[3];
	idx_t *xadj, *vwgt, *adjncy, *adjwgt, *label, *where, *bndptr, *bndind;
	idx_t *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *slabel[2];
	idx_t *rename;
	idx_t *auxadjncy;
	graph_t *lgraph, *rgraph, *lrborder, *rrborder, *lcborder, *rcborder, *corner;

	WCOREPUSH;

	IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->SplitTmr));

	nvtxs   = graph->nvtxs;
	xadj    = graph->xadj;
	vwgt    = graph->vwgt;
	adjncy  = graph->adjncy;
	adjwgt  = graph->adjwgt;
	label   = graph->label;
	where   = graph->where;
	bndptr  = graph->bndptr;
	bndind  = graph->bndind;
	ASSERT(bndptr != NULL);

	rename = iwspacemalloc(ctrl, nvtxs);

	snvtxs[0] = snvtxs[1] = snvtxs[2] = snedges[0] = snedges[1] = snedges[2] = 0;
	for (i = 0; i < nvtxs; i++) {
		k = where[i];
		rename[i] = snvtxs[k]++;	/* The vertex index of each vertex in a subgraph */
		snedges[k] += xadj[i+1]-xadj[i];	/* XXX should not removed edges be deleted ? */
	}

	lgraph      = SetupSplitGraph(graph, snvtxs[0], snedges[0]);
	sxadj[0]    = lgraph->xadj;
	svwgt[0]    = lgraph->vwgt;
	sadjncy[0]  = lgraph->adjncy;
	sadjwgt[0]  = lgraph->adjwgt;
	slabel[0]   = lgraph->label;

	rgraph      = SetupSplitGraph(graph, snvtxs[1], snedges[1]);
	sxadj[1]    = rgraph->xadj;
	svwgt[1]    = rgraph->vwgt;
	sadjncy[1]  = rgraph->adjncy;
	sadjwgt[1]  = rgraph->adjwgt;
	slabel[1]   = rgraph->label;

	/* Go and use bndptr to also mark the boundary nodes in the two partitions */
	for (ii = 0; ii < graph->nbnd; ii++) {
		i = bndind[ii];
		for (j = xadj[i]; j < xadj[i+1]; j++)
			bndptr[adjncy[j]] = 1;
	}

	/* restat snvtxs and snedges */
	snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
	sxadj[0][0] = sxadj[1][0] = 0;
	for (i = 0; i < nvtxs; i++) {
		if ((mypart = where[i]) == 2)
			continue;

		istart = xadj[i];
		iend   = xadj[i+1];
		if (bndptr[i] == -1) { /* This is an interior vertex */
			/* auxadjncy should be the start position to write adjacencies,
			 * the reason for minus 'istart' is that the for loop starts from 'istart' */
			auxadjncy = sadjncy[mypart] + snedges[mypart] - istart;
			for(j = istart; j < iend; j++)
				auxadjncy[j] = adjncy[j];
			snedges[mypart] += iend - istart;
		}
		else {	/* This is a boundary vertex */
			auxadjncy = sadjncy[mypart];
			l = snedges[mypart];
			for (j = istart; j < iend; j++) {
				k = adjncy[j];
				if (where[k] == mypart)
					auxadjncy[l++] = k;
			}
			snedges[mypart] = l;
		}

		svwgt[mypart][snvtxs[mypart]]    = vwgt[i];
		slabel[mypart][snvtxs[mypart]]   = label[i];
		sxadj[mypart][++snvtxs[mypart]]  = snedges[mypart];
	}

	for (mypart = 0; mypart < 2; mypart++) {
		iend = snedges[mypart];
		iset(iend, 1, sadjwgt[mypart]);

		auxadjncy = sadjncy[mypart];
		for (i = 0; i < iend; i++)
			auxadjncy[i] = rename[auxadjncy[i]];
	}

	lgraph->nvtxs  = snvtxs[0];
	lgraph->nedges = snedges[0];
	rgraph->nvtxs  = snvtxs[1];
	rgraph->nedges = snedges[1];

	SetupGraph_tvwgt(lgraph);
	SetupGraph_tvwgt(rgraph);

	IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->SplitTmr));

	*r_lgraph = lgraph;
	*r_rgraph = rgraph;

	WCOREPOP;
}

/**
 * This function constructs lgraph and rgraph for a compressed graph 'before'
 * according to the separation results of 'after' and the relationship between
 * 'before' and 'after'. Note that the final lgraph and rgraph constructed
 * correspond to the original graph 'before'
 * XXX new nbnd, bndind, label of sort[i]->super should be modified in this function.
 */
void SplitGraphOrderUncompressBDF(ctrl_t *ctrl, graph_t *graph, graph_t *cgraph, idx_t *cptr,
		idx_t *cind, graph_t **r_lgraph, graph_t **r_rgraph)
{
	/* construct perm from iperm */
	/*
	for (i=0; i<nnvtxs; i++)
		perm[iperm[i]] = i;
	for (l=ii=0; ii<nnvtxs; ii++) {
		i = perm[ii];
		for (j=cptr[i]; j<cptr[i+1]; j++)
			iperm[cind[j]] = l++;
	}
	gk_free((void **)&cptr, &cind, LTERM);
	*/
}

/**
 * This function extracts a bigraph from the original graph according
 * to a subgraph (lgraph/rgraph) and boundary nodes.
 * Note: # of boundary nodes is nbnd = graph->nbnd
 *		indices of boundary is bndind = graph->bndind
 *		node label is stored in label = graph->label
 *	Namely, all necessary partitioning information has been prepared in graph.
 * 	/param graph is the original graph
 * 	/param sgraph is the subgraph split from the original graph
 * 	/param nrows is the total number of rows in the WHOLE original graph
 */
void ExtractBiGraph(ctrl_t *ctrl, bigraph_t *bigraph, graph_t *lgraph, graph_t *rgraph,
		bigraph_t **r_lbigraph, bigraph_t **r_rbigraph)
{
	bigraph_t *lbigraph = NULL, *rbigraph = NULL, *border = NULL;
	bigraph_t *p, *q, *r, *s;
	graph_t *graph = bigraph->super;
	idx_t i, rj, cj;
	idx_t lnrows = 0, lncols = 0;
	idx_t rnrows = 0, rncols = 0;
	idx_t nrbnds = 0, ncbnds = 0;
	label_t *lrlabel, *lclabel, *rrlabel, *rclabel, *rblabel, *cblabel;
	idx_t debugnz = 0, nz, area, snz, sarea;

	/* stat # nodes in each block */
	for (i = 0; i < lgraph->nvtxs; i++) {
		if (lgraph->label[i] < ctrl->nrows)	lnrows++;
		else	lncols++;
	}
	for (i = 0; i < rgraph->nvtxs; i++) {
		if (rgraph->label[i] < ctrl->nrows)	rnrows++;
		else	rncols++;
	}
	for (i = 0; i < graph->nbnd; i++) {
		if (graph->label[graph->bndind[i]] < ctrl->nrows)	nrbnds++;
		else	ncbnds++;
	}
	/* allocate memory for commonly used labels
	 * NOTE ! We use the same label for all borders in order to save memory! */
	lrlabel = (label_t*)gk_malloc(sizeof(label_t), "ExtractBiGraph : lrlabel");
	lclabel = (label_t*)gk_malloc(sizeof(label_t), "ExtractBiGraph : lclabel");
	rrlabel = (label_t*)gk_malloc(sizeof(label_t), "ExtractBiGraph : rrlabel");
	rclabel = (label_t*)gk_malloc(sizeof(label_t), "ExtractBiGraph : rclabel");
	rblabel = (label_t*)gk_malloc(sizeof(label_t), "ExtractBiGraph : rblabel");
	cblabel = (label_t*)gk_malloc(sizeof(label_t), "ExtractBiGraph : cblabel");

	lrlabel->label = imalloc(lnrows, "ExtractBiGraph : lrlabel->label");
	lclabel->label = imalloc(lncols, "ExtractBiGraph : lclabel->label");
	rrlabel->label = imalloc(rnrows, "ExtractBiGraph : rrlabel->label");
	rclabel->label = imalloc(rncols, "ExtractBiGraph : rclabel->label");
	rblabel->label = imalloc(nrbnds, "ExtractBiGraph : rblabel->label");
	cblabel->label = imalloc(ncbnds, "ExtractBiGraph : cblabel->label");

	lrlabel->ref = lclabel->ref = rrlabel->ref = rclabel->ref = rblabel->ref = cblabel->ref = 0;

	for (rj = cj = i = 0; i < lgraph->nvtxs; i++) {
		if (lgraph->label[i] < ctrl->nrows)	lrlabel->label[rj++] = lgraph->label[i];
		else	lclabel->label[cj++] = lgraph->label[i];
	}
	for (rj = cj = i = 0; i < rgraph->nvtxs; i++) {
		if (rgraph->label[i] < ctrl->nrows)	rrlabel->label[rj++] = rgraph->label[i];
		else	rclabel->label[cj++] = rgraph->label[i];
	}
	for (rj = cj = i = 0; i < graph->nbnd; i++) {
		if (graph->label[graph->bndind[i]] < ctrl->nrows)	rblabel->label[rj++] = graph->label[graph->bndind[i]];
		else	cblabel->label[cj++] = graph->label[graph->bndind[i]];
	}

	gk_startcputimer(_nztimer);	/* TODO exp timer */

	/* construct two bigraphs */
	lbigraph = SetupBiGraphFromGraph(lgraph, lnrows, lncols, ctrl->nrows, ctrl->ncols, lrlabel, lclabel);

	/* construct the fist column border */
	border = CreateBorder(lnrows, ncbnds, StatNonZeros(ctrl, lrlabel, cblabel, lnrows, ncbnds), lrlabel, cblabel);
	lbigraph->right = border;

	/* split original column borders for new column borders */
	p = border;
	q = bigraph->right;
	while (q) {
		border = CreateBorder(lnrows, q->ncols, StatNonZeros(ctrl, lrlabel, q->clabel, lnrows, q->ncols), lrlabel, q->clabel);
		p->right = border;
		p = p->right;
		q = q->right;
	}

	/* construct the first row border */
	border = CreateBorder(nrbnds, lncols, StatNonZeros(ctrl, rblabel, lclabel, nrbnds, lncols), rblabel, lclabel);
	lbigraph->down = border;

	/* construct the conner */
	p = border;
	border = CreateBorder(nrbnds, ncbnds, StatNonZeros(ctrl, rblabel, cblabel, nrbnds, ncbnds), rblabel, cblabel);
	debugnz = border->nz;
	p->right = border;

	/* split original column borders for new conners */
	p = border;
	q = bigraph->right;
	while (q) {
		border = CreateBorder(nrbnds, q->ncols, StatNonZeros(ctrl, rblabel, q->clabel, nrbnds, q->ncols), rblabel, q->clabel);
		p->right = border;
		p = p->right;
		q = q->right;
	}

	/* split original row borders for new row borders/corners,
	 * connect orginal cornners to the new bigraph */
	p = lbigraph->down;
	q = bigraph->down;
	while (q) {
		/* row border from original row border */
		border = CreateBorder(q->nrows, lncols, StatNonZeros(ctrl, q->rlabel, lclabel, q->nrows, lncols), q->rlabel, lclabel);
		p->down = border;
		p = p->down;
		/* row conner from original row border */
		p->right = CreateBorder(q->nrows, ncbnds, StatNonZeros(ctrl, q->rlabel, cblabel, q->nrows, ncbnds), q->rlabel, cblabel);
		/* connect the original cornners */
		s = p->right;
		r = q->right;
		while (r) {
			s->right = CreateBorder(r->nrows, r->ncols, r->nz, r->rlabel, r->clabel);
			r->rlabel->ref++;	r->clabel->ref++;
			s = s->right;
			r = r->right;
		}
		q = q->down;
	}

	/************************************************************/
	/* construct the right bigraph */
	rbigraph = SetupBiGraphFromGraph(rgraph, rnrows, rncols, ctrl->nrows, ctrl->ncols, rrlabel, rclabel);

	/* construct the fist column border */
	border = CreateBorder(rnrows, ncbnds, StatNonZeros(ctrl, rrlabel, cblabel, rnrows, ncbnds), rrlabel, cblabel);
	rbigraph->right = border;

	/* split original column borders for new column borders */
	p = border;
	q = bigraph->right;
	while (q) {
		border = CreateBorder(rnrows, q->ncols, StatNonZeros(ctrl, rrlabel, q->clabel, rnrows, q->ncols), rrlabel, q->clabel);
		p->right = border;
		p = p->right;
		q = q->right;
	}

	/* construct the first row border */
	border = CreateBorder(nrbnds, rncols, StatNonZeros(ctrl, rblabel, rclabel, nrbnds, rncols), rblabel, rclabel);
	rbigraph->down = border;

	/* construct the conner */
	p = border;
	border = CreateBorder(nrbnds, ncbnds, StatNonZeros(ctrl, rblabel, cblabel, nrbnds, ncbnds), rblabel, cblabel);
	p->right = border;

	/* split original column borders for new conners */
	p = border;
	q = bigraph->right;
	while (q) {
		border = CreateBorder(nrbnds, q->ncols, StatNonZeros(ctrl, rblabel, q->clabel, nrbnds, q->ncols), rblabel, q->clabel);
		p->right = border;
		p = p->right;
		q = q->right;
	}

	/* split original row borders for new row borders/corners,
	 * connect orginal cornners to the new bigraph */
	p = rbigraph->down;
	q = bigraph->down;
	while (q) {
		/* row border from original row border */
		border = CreateBorder(q->nrows, rncols, StatNonZeros(ctrl, q->rlabel, rclabel, q->nrows, rncols), q->rlabel, rclabel);
		p->down = border;
		p = p->down;
		/* row conner from original row border */
		p->right = CreateBorder(q->nrows, ncbnds, StatNonZeros(ctrl, q->rlabel, cblabel, q->nrows, ncbnds), q->rlabel, cblabel);
		/* connect the original cornners */
		s = p->right;
		r = q->right;
		while (r) {
			s->right = CreateBorder(r->nrows, r->ncols, r->nz, r->rlabel, r->clabel);
			s = s->right;
			r = r->right;
		}
		q = q->down;
	}

	gk_stopcputimer(_nztimer);	/* TODO exp timer */

	*r_lbigraph = lbigraph;
	*r_rbigraph = rbigraph;

	ASSERT(CheckArea(bigraph, lbigraph, rbigraph));	/* TODO debug */
	ASSERT(CheckNonZeros(bigraph, lbigraph, rbigraph));	/*TODO debug*/
}

void swap(idx_t *x,idx_t *y)
{
   idx_t temp;
   temp = *x;
   *x = *y;
   *y = temp;
}

int choose_pivot(int i,int j )
{
   return((i+j) /2);
}

void quicksort(idx_t list[],int m,int n)
{
   int i,j,k;
   idx_t key;
   if( m < n)
   {
      k = choose_pivot(m,n);
      swap(&list[m],&list[k]);
      key = list[m];
      i = m+1;
      j = n;
      while(i <= j)
      {
         while((i <= n) && (list[i] <= key))
                i++;
         while((j >= m) && (list[j] > key))
                j--;
         if( i < j)
                swap(&list[i],&list[j]);
      }
      swap(&list[m],&list[j]);
      quicksort(list,m,j-1);
      quicksort(list,j+1,n);
   }
}

idx_t StatNonZeros(ctrl_t *ctrl, label_t *rlabel, label_t *clabel, idx_t nrows, idx_t ncols) {
	idx_t i, j, k;
	idx_t istart, iend, adj;
	idx_t nz = 0;
	if (nrows < ncols) {
		idx_t *tmp = (idx_t*)malloc(sizeof(idx_t)*ncols);
	        for (i = 0; i < ncols; i++)
				tmp[i] = clabel->label[i];
		quicksort(tmp, 0, ncols-1);	
		for (i = 0; i < nrows; i++) {
			istart = ctrl->obigraph->super->xadj[rlabel->label[i]];
			iend = ctrl->obigraph->super->xadj[rlabel->label[i]+1];
			for (j = istart; j < iend; j++) {
				adj = ctrl->obigraph->super->adjncy[j];
				int head = 0, tail = ncols, mid = (head + tail) / 2;
				while (head < tail-1)
				{
					if (tmp[mid] <= adj) head = mid;
					else tail = mid;
					mid = (head + tail) / 2;
				}
				if (tmp[head] == adj)
					nz++;
			}
		}
		free(tmp);
	}
	else {
		idx_t *tmp = (idx_t*)malloc(sizeof(idx_t)*nrows);
	        for (i = 0; i < nrows; i++)
			tmp[i] = rlabel->label[i];
		quicksort(tmp, 0, nrows-1);	
		for (i = 0; i < ncols; i++) {
			istart = ctrl->obigraph->super->xadj[clabel->label[i]];
			iend = ctrl->obigraph->super->xadj[clabel->label[i]+1];
			for (j = istart; j < iend; j++) {
				adj = ctrl->obigraph->super->adjncy[j];
				int head = 0, tail = nrows, mid = (head + tail) / 2;
				int o1=0, o2= 0;
				while (head < tail-1)
				{
					if (tmp[mid] <= adj) head = mid;
					else tail = mid;
					mid = (head + tail) / 2;
				}
				if (tmp[head] == adj)
				{
					nz++;
				}
			}
		}
		free(tmp);
	}

	return nz;
}

/**
 * This function checks whether the area after partitioning is equal
 * to that before partitioning.
 */
idx_t CheckArea(bigraph_t *bigraph, bigraph_t *lbigraph, bigraph_t *rbigraph){
	idx_t sarea = 0, sareabytimes = 0;

	sarea += lbigraph->area + lbigraph->right->area + lbigraph->down->area + lbigraph->down->right->area;
	sarea += rbigraph->area + rbigraph->right->area + rbigraph->down->area + rbigraph->down->right->area;
	sarea -= lbigraph->down->right->area;
	sarea += lbigraph->nrows * rbigraph->ncols + rbigraph->nrows * lbigraph->ncols;

	if (sarea != bigraph->area)	return 0;
	if (sarea != bigraph->nrows * bigraph->ncols) return 0;

	return 1;
}

/**
 * This function checks whether the total non-zeros after partitioning is equal
 * to that before partitioning.
 */
int CheckNonZeros(bigraph_t *bigraph, bigraph_t *lbigraph, bigraph_t *rbigraph) {
	idx_t onz = 0, nnz = 0, nz;
	idx_t area;
	bigraph_t *p, *q;

	StatNzAndArea(bigraph, &nz, &area, 0);
	onz = nz;

	StatNzAndArea(lbigraph, &nz, &area, 0);
	nnz += nz;
	StatNzAndArea(rbigraph, &nz, &area, 0);
	nnz += nz;

	nz = 0;
	p = lbigraph->down;
	while (p) {
		q = p->right;
		while (q) {
			nz += q->nz;
			q = q->right;
		}
		p = p->down;
	}
	nnz -= nz;

	if (nnz != onz) {
		printf("***ERROR : non-zero error !\n");
		return 0;
	}
	return 1;
}

/**
 * This function checks the format of final perm and iperm.
 */
idx_t CheckPermIPerm(idx_t *perm, idx_t *iperm, idx_t nvtxs) {
	idx_t i, psum = 0, ipsum = 0, realsum = 0;
	for (i = 0; i < nvtxs; i++) {
		psum += perm[i];
		ipsum += iperm[i];
	}
	realsum = (nvtxs - 1) * nvtxs / 2;	/* from 0 on */
	if (ipsum != realsum)	return 0;
	if (psum != realsum)	return 0;
	return 1;
}

/**
 * This function prints the sorted graph list by their area of density.
 */
void PrintSortedList(idx_t ndiags, bigraph_t **sort, idx_t *areas, real_t *denses, idx_t isarea) {
	idx_t i;

	if (isarea) {
		printf("--tareas list: ");
		for (i = 0; i < ndiags; i++)	printf("%lld;\t", areas[i]);
		printf("\n");
	}
	else {
		printf("---denses list: ");
		for (i = 0; i < ndiags; i++)	printf("%.6f;\t", denses[i]);
		printf("\n");
	}

	printf("---areas list: ");
	for (i = 0; i < ndiags; i++)	printf("%lld;\t", sort[i]->area);
	printf("\n---nrows list: ");
	for (i = 0; i < ndiags; i++)	printf("%d;\t", sort[i]->nrows);
	printf("\n---ncols list: ");
	for (i = 0; i < ndiags; i++)	printf("%d;\t", sort[i]->ncols);
	printf("\n");
}
