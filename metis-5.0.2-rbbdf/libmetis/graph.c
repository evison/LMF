/**
\file
\brief Functions that deal with setting up the graphs for METIS.

\date   Started 7/25/1997
\author George  
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version\verbatim $Id: graph.c 10513 2011-07-07 22:06:03Z karypis $ \endverbatim
*/

#include "metislib.h"


/*************************************************************************/
/*! This function sets up the graph from the user input */
/*************************************************************************/
graph_t *SetupGraph(ctrl_t *ctrl, idx_t nvtxs, idx_t ncon, idx_t *xadj, 
             idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt) 
{
  idx_t i, j, k, sum;
  real_t *nvwgt;
  graph_t *graph;

  /* allocate the graph and fill in the fields */
  graph = CreateGraph();

  graph->nvtxs  = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon   = ncon;

  graph->xadj      = xadj;
  graph->free_xadj = 0;

  graph->adjncy      = adjncy;
  graph->free_adjncy = 0;


  /* setup the vertex weights */
  if (vwgt) {
    graph->vwgt      = vwgt;
    graph->free_vwgt = 0;
  }
  else {
    vwgt = graph->vwgt = ismalloc(ncon*nvtxs, 1, "SetupGraph: vwgt");
  }

  graph->tvwgt    = imalloc(ncon, "SetupGraph: tvwgts");
  graph->invtvwgt = rmalloc(ncon, "SetupGraph: invtvwgts");
  for (i=0; i<ncon; i++) {
    graph->tvwgt[i]    = isum(nvtxs, vwgt+i, ncon);
    graph->invtvwgt[i] = 1.0/(graph->tvwgt[i] > 0 ? graph->tvwgt[i] : 1);
  }


  if (ctrl->objtype == METIS_OBJTYPE_VOL) { 
    /* Setup the vsize */
    if (vsize) {
      graph->vsize      = vsize;
      graph->free_vsize = 0;
    }
    else {
      vsize = graph->vsize = ismalloc(nvtxs, 1, "SetupGraph: vsize");
    }

    /* Allocate memory for edge weights and initialize them to the sum of the vsize */
    adjwgt = graph->adjwgt = imalloc(graph->nedges, "SetupGraph: adjwgt");
    for (i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        adjwgt[j] = 1+vsize[i]+vsize[adjncy[j]];
    }
  }
  else { /* For edgecut minimization */
    /* setup the edge weights */
    if (adjwgt) {
      graph->adjwgt      = adjwgt;
      graph->free_adjwgt = 0;
    }
    else {
      adjwgt = graph->adjwgt = ismalloc(graph->nedges, 1, "SetupGraph: adjwgt");
    }
  }


  /* setup various derived info */
  SetupGraph_tvwgt(graph);

  /* evison */
  if (ctrl->optype == METIS_OP_PMETIS || ctrl->optype == METIS_OP_OMETIS || ctrl->optype == METIS_OP_BMETIS)
    SetupGraph_label(graph);

  ASSERT(CheckGraph(graph, ctrl->numflag, 1));

  return graph;
}

/* evison */
bigraph_t* SetupBiGraphFromParams(ctrl_t *ctrl, idx_t nvtxs, idx_t ncon, idx_t *xadj,
             idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
             idx_t nrows, idx_t ncols, idx_t totalnrows, idx_t totalncols,
             label_t *rlabel, label_t *clabel){
	graph_t *graph;
	bigraph_t *bigraph;

	if (nvtxs != nrows + ncols)	return NULL;

	graph = SetupGraph(ctrl, nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt);

	if (!graph)	return NULL;

	bigraph = SetupBiGraphFromGraph(graph, nrows, ncols, totalnrows,
			totalncols, rlabel, clabel);

	return bigraph;
}

bigraph_t* SetupBiGraphFromGraph(graph_t *graph, idx_t nrows, idx_t ncols,
		idx_t totalnrows, idx_t totalncols, label_t *rlabel, label_t *clabel) {
	bigraph_t *bigraph;
	int i, j, k;

	if(!graph)	return NULL;
	if(graph->nvtxs != nrows + ncols){
		printf("***ERROR bigraph detected !\n");
		return NULL;
	}

	bigraph = CreateBiGraph();

	bigraph->super = graph;
	bigraph->lastvtx = graph->nvtxs;
	bigraph->nrows = nrows;
	bigraph->ncols = ncols;
	bigraph->area = nrows * ncols;
	bigraph->nz = graph->nedges / 2;
	bigraph->partible = 1;

	bigraph->rlabel = rlabel;
	bigraph->clabel = clabel;
	bigraph->rlabel->ref++;
	bigraph->clabel->ref++;
	/*
	bigraph->rlabel = imalloc(bigraph->nrows, "SetupBiGraph_rlabel: rlabel");
	bigraph->clabel = imalloc(bigraph->ncols, "SetupBiGraph_clabel: clabel");
	i = j = 0;
	for (k = 0; k < graph->nvtxs; k++) {
		if (graph->label[k] < totalnrows)	bigraph->rlabel[i++] = graph->label[k];
		else	bigraph->clabel[j++] = graph->label[k];
	}
	*/
	/*
	bigraph->rlabel = (label_t*)gk_malloc(sizeof(label_t), "SetupBiGraphFromGraph : bigraph->rlabel");
	bigraph->rlabel->label = imalloc(bigraph->nrows, "SetupBiGraphFromGraph: bigraph->rlabel->label");
	bigraph->clabel = (label_t*)gk_malloc(sizeof(label_t), "SetupBiGraphFromGraph : bigraph->clabel");
	bigraph->clabel->label = imalloc(bigraph->ncols, "SetupBiGraphFromGraph: bigraph->clabel->label");
	i = j = 0;
	for (k = 0; k < graph->nvtxs; k++) {
		if (graph->label[k] < totalnrows)	bigraph->rlabel->label[i++] = graph->label[k];
		else	bigraph->clabel->label[j++] = graph->label[k];
	}
	bigraph->rlabel->ref = 1;
	bigraph->clabel->ref = 1;
	*/

	bigraph->right = NULL;
	bigraph->down = NULL;
	bigraph->next = NULL;

	return bigraph;
}


/*************************************************************************/
/*! Set's up the tvwgt/invtvwgt info */
/*************************************************************************/
void SetupGraph_tvwgt(graph_t *graph)
{
  idx_t i;

  if (graph->tvwgt == NULL) 
    graph->tvwgt  = imalloc(graph->ncon, "SetupGraph_tvwgt: tvwgt");
  if (graph->invtvwgt == NULL) 
    graph->invtvwgt = rmalloc(graph->ncon, "SetupGraph_tvwgt: invtvwgt");

  for (i=0; i<graph->ncon; i++) {
    graph->tvwgt[i]    = isum(graph->nvtxs, graph->vwgt+i, graph->ncon);
    graph->invtvwgt[i] = 1.0/(graph->tvwgt[i] > 0 ? graph->tvwgt[i] : 1);
  }
}


/*************************************************************************/
/*! Set's up the label info */
/*************************************************************************/
void SetupGraph_label(graph_t *graph)
{
  idx_t i;

  if (graph->label == NULL)
    graph->label = imalloc(graph->nvtxs, "SetupGraph_label: label");

  for (i=0; i<graph->nvtxs; i++)
    graph->label[i] = i;	/* TODO from 0 or 1 on ? */
}


/*************************************************************************/
/*! Setup the various arrays for the splitted graph */
/*************************************************************************/
graph_t *SetupSplitGraph(graph_t *graph, idx_t snvtxs, idx_t snedges)
{
  graph_t *sgraph;

  sgraph = CreateGraph();

  sgraph->nvtxs  = snvtxs;
  sgraph->nedges = snedges;
  sgraph->ncon   = graph->ncon;

  /* Allocate memory for the splitted graph */
  sgraph->xadj        = imalloc(snvtxs+1, "SetupSplitGraph: xadj");
  sgraph->vwgt        = imalloc(sgraph->ncon*snvtxs, "SetupSplitGraph: vwgt");
  sgraph->adjncy      = imalloc(snedges,  "SetupSplitGraph: adjncy");
  sgraph->adjwgt      = imalloc(snedges,  "SetupSplitGraph: adjwgt");
  sgraph->label	      = imalloc(snvtxs,   "SetupSplitGraph: label");
  sgraph->tvwgt       = imalloc(sgraph->ncon, "SetupSplitGraph: tvwgt");
  sgraph->invtvwgt    = rmalloc(sgraph->ncon, "SetupSplitGraph: invtvwgt");

  if (graph->vsize)
    sgraph->vsize     = imalloc(snvtxs,   "SetupSplitGraph: vsize");

  return sgraph;
}


/*************************************************************************/
/*! This function creates and initializes a graph_t data structure */
/*************************************************************************/
graph_t *CreateGraph(void)
{
  graph_t *graph;

  graph = (graph_t *)gk_malloc(sizeof(graph_t), "CreateGraph: graph");

  InitGraph(graph);

  return graph;
}

/*evison*/
bigraph_t *CreateBiGraph(void){
	bigraph_t *bigraph;

	bigraph = (bigraph_t *)gk_malloc(sizeof(bigraph_t), "CreateBiGraph: bigraph");

	InitBiGraph(bigraph);

	return bigraph;
}

bigraph_t *CreateBorder(idx_t nrows, idx_t ncols, idx_t nz, label_t *rlabel, label_t *clabel) {
	bigraph_t *border;

	border = CreateBiGraph();
	border->nrows = nrows;
	border->ncols = ncols;
	border->nz = nz;
	border->area = nrows * ncols;
	border->rlabel = rlabel;
	border->clabel = clabel;
	border->rlabel->ref++;
	border->clabel->ref++;
	border->partible = 1;

	return border;
}

/*************************************************************************/
/*! This function initializes a graph_t data structure */
/*************************************************************************/
void InitGraph(graph_t *graph) 
{
  memset((void *)graph, 0, sizeof(graph_t));

  /* graph size constants */
  graph->nvtxs     = -1;
  graph->nedges    = -1;
  graph->ncon      = -1;
  graph->mincut    = -1;
  graph->minvol    = -1;
  graph->nbnd      = -1;

  /* memory for the graph structure */
  graph->xadj      = NULL;
  graph->vwgt      = NULL;
  graph->vsize     = NULL;
  graph->adjncy    = NULL;
  graph->adjwgt    = NULL;
  graph->label     = NULL;
  graph->cmap      = NULL;
  graph->tvwgt     = NULL;
  graph->invtvwgt  = NULL;

  /* by default these are set to true, but the can be explicitly changed afterwards */
  graph->free_xadj   = 1;
  graph->free_vwgt   = 1;
  graph->free_vsize  = 1;
  graph->free_adjncy = 1;
  graph->free_adjwgt = 1;


  /* memory for the partition/refinement structure */
  graph->where     = NULL;
  graph->pwgts     = NULL;
  graph->id        = NULL;
  graph->ed        = NULL;
  graph->bndptr    = NULL;
  graph->bndind    = NULL;
  graph->nrinfo    = NULL;
  graph->ckrinfo   = NULL;
  graph->vkrinfo   = NULL;

  /* linked-list structure */
  graph->coarser   = NULL;
  graph->finer     = NULL;
}

void InitBiGraph(bigraph_t *bigraph){
	memset((void *)bigraph, 0, sizeof(bigraph_t));
	bigraph->super = NULL;
	bigraph->lastvtx = -1;
	bigraph->ncols = -1;
	bigraph->nrows = -1;
	bigraph->area = -1;
	bigraph->nz = -1;
	bigraph->partible = 1;
	bigraph->rlabel = NULL;
	bigraph->clabel= NULL;
	bigraph->right = NULL;
	bigraph->down = NULL;
	bigraph->next = NULL;
}

/*************************************************************************/
/*! This function frees the refinement/partition memory stored in a graph */
/*************************************************************************/
void FreeRData(graph_t *graph) 
{

  /* The following is for the -minconn and -contig to work properly in
     the vol-refinement routines */
  if ((void *)graph->ckrinfo == (void *)graph->vkrinfo)
    graph->ckrinfo = NULL;

  /* free partition/refinement structure */
  gk_free((void **)&graph->where, &graph->pwgts, &graph->id, &graph->ed, 
      &graph->bndptr, &graph->bndind, &graph->nrinfo, &graph->ckrinfo, 
      &graph->vkrinfo, LTERM);
}


/*************************************************************************/
/*! This function deallocates any memory stored in a graph */
/*************************************************************************/
void FreeGraph(graph_t **r_graph) 
{
  graph_t *graph;

  graph = *r_graph;

  /* free graph structure */
  if (graph->free_xadj){
    gk_free((void **)&graph->xadj, LTERM);
  }

  if (graph->free_vwgt)
    gk_free((void **)&graph->vwgt, LTERM);

  if (graph->free_vsize)
    gk_free((void **)&graph->vsize, LTERM);

  if (graph->free_adjncy)
    gk_free((void **)&graph->adjncy, LTERM);

  if (graph->free_adjwgt)
    gk_free((void **)&graph->adjwgt, LTERM);

  /* free partition/refinement structure */
  FreeRData(graph);

  gk_free((void **)&graph->tvwgt, &graph->invtvwgt, &graph->label, 
      &graph->cmap, &graph, LTERM);

  *r_graph = NULL;
}


/*******************************************************************************
 * This function frees a bipartite graph.
 * This function CAN ONLY be used by FreeBigraphBorder.
 *******************************************************************************/
void FreeBiGraph(ctrl_t *ctrl, bigraph_t **r_bigraph){
	bigraph_t *bigraph = *r_bigraph;

	if (ctrl == NULL || bigraph == NULL)
		return;

	if (bigraph == ctrl->obigraph)
		return;

	if (bigraph->super)	FreeGraph(&bigraph->super);

	if (bigraph->rlabel){
		bigraph->rlabel->ref--;
		if (bigraph->rlabel->ref == 0)	gk_free((void**)&bigraph->rlabel->label, LTERM);
		bigraph->rlabel = NULL;
	}
	if (bigraph->clabel){
		bigraph->clabel->ref--;
		if (bigraph->clabel->ref == 0)	gk_free((void**)&bigraph->clabel->label, LTERM);
		bigraph->clabel = NULL;
	}

	/* Note: we just ignore the free of pointers right, down and next here */

	gk_free((void**)r_bigraph, LTERM);	/* TODO */
	*r_bigraph = NULL;
}

/*******************************************************************************
 * This function frees a block diagonal as well as the borders realated to it.
 * This function can be used outside.
 *******************************************************************************/
void FreeBiGraphBorder(ctrl_t *ctrl, bigraph_t **r_biborder) {
	bigraph_t *p, *q, *r = *r_biborder;

	while (r) {
		p = r;
		while (p) {
			q = p->right;
			FreeBiGraph(ctrl, &p);
			p = q;
		}
		r = r->down;
	}
	*r_biborder = NULL;
}

void FreeBiGraphBorderList(ctrl_t *ctrl, bigraph_t **r_biborderlist, idx_t ndiags){
	bigraph_t *p = *r_biborderlist, *q;
	idx_t i = 0;

	if(p == NULL)	return;

	printf("%d diagonal block to free.\n", ndiags);

	while(p){
		i++;
		printf("freeing D%d\r", i);
		q = p->next;
		FreeBiGraphBorder(ctrl, &p);
		p = q;
	}
	printf("\n");

	*r_biborderlist = NULL;
}
