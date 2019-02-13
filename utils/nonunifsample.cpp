#include "mex.h"
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <iostream>

static double* prefix_sums = NULL;
static mwSize* num_elem = NULL;
static mwSize* tree_size = NULL;

void exitFcn()
{
	if( num_elem != NULL )
	  mxFree(num_elem);

	if( prefix_sums != NULL )
	  mxFree(prefix_sums);

	if( tree_size != NULL )
      mxFree(tree_size);
}

void build_prefix_tree(double* weights, int num_elem)
{
  for( int k=num_elem-1; k >= 0; k-- )
  {
  	prefix_sums[(*tree_size/2)+k] = std::max(weights[k],0.0);
  }

  for( int k=num_elem-1; k >= 0; k-- )
  {
  	prefix_sums[k] = prefix_sums[2*k+1] + prefix_sums[2*k+2];
  }
}

int sample_helper(double r, int k)
{
	// Leaf of tree
	if( 2*k+1 >= *tree_size )
	  return k;

	double left = prefix_sums[2*k+1];

	if( r <= left )
	  return sample_helper(r, 2*k+1);
	else
	  return sample_helper(r - left, 2*k+2);
}

int sample(double r)
{
  // Rescale uniformly random variable
  r *= prefix_sums[0];

  int k=0;
  return sample_helper(r,k) - (*tree_size)/2+1;
}

void update_weight( double w, int k )
{
  int j = (*tree_size)/2 + k - 1; // Shift by 1 for 0-based indexing
  prefix_sums[j] = std::max(w,0.0);

  j = (j-1)/2;
  while( j >= 0 )
  {
  	prefix_sums[j] = prefix_sums[2*j+1] + prefix_sums[2*j+2];
  	if( j==0 )
  		break;
  	j = (j-1)/2;
  }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  // Get the command string.
  char cmd[64];
  if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
    mexErrMsgTxt("First input should be a command string less than 64 characters long.");

  // Command: Init
  if( !strcmp("init", cmd) )
  {
	mexAtExit(exitFcn);

    num_elem = (mwSize*) mxCalloc(1,sizeof(mwSize));
    mexMakeMemoryPersistent(num_elem);

    tree_size = (mwSize*) mxCalloc(1,sizeof(mwSize));
    mexMakeMemoryPersistent(tree_size);

    mwSize m = mxGetM(prhs[1]);
    *num_elem = m;

    *tree_size = (2 << uint(ceil(log2(double(m)))))-1;

    prefix_sums = (double*) mxCalloc(*tree_size, sizeof(double));
    mexMakeMemoryPersistent(prefix_sums);

    double* weights = mxGetPr(prhs[1]);

    build_prefix_tree(weights, *num_elem);
  }

  if( !strcmp("sample", cmd) )
  {
  	double r = mxGetScalar(prhs[1]);
  	int k = sample(r);

  	plhs[0] = mxCreateDoubleScalar(k);
  }

  if( !strcmp("update", cmd) )
  {
  	double w = mxGetScalar(prhs[1]);
  	int k = int(mxGetScalar(prhs[2]));
  	update_weight(w, k);
  }

  return;
}