#include <math.h>

/************************************************************************
 * Function int _get2dPhase                                             *
 * Computes the integrated phase along one given direction (angle)      *
 * All init data are pre-computed to accelerate execution time within   *
 * the time critical aoloop() in yao.i.                                 *
 * Last modified: Dec 12, 2003.                                         *
 * Author: F.Rigaut                                                     *
 ************************************************************************/

int _dist(float *d, long dimx, long dimy, float xc, float yc)
{
  /* Declarations */
  long i,j;

  /* Loop and fill d with distance values */
  for (i=0;i<dimx;++i) {
    for (j=0;j<dimy;++j) {
      d[i + j * dimx] = (float)sqrt( (xc-(float)i) * (xc-(float)i) + 
			 	     (yc-(float)j) * (yc-(float)j) );
    }
  }
  return 0;
}

int _get2dPhase(float *pscreens, /* dimension [psnx,psny,nscreens] */
    int psnx, 
    int psny, 
    int nscreens, 
    float *outphase, /* dimension [phnx,phny] */
    int phnx, 
    int phny, 
    int *ishifts,    /* array of X integer shifts dimension [phnx,nscreens] */
    float *xshifts,  /* array of X fractional shifts dimension [phnx,nscreens] */
    int *jshifts,    /* array of Y integer shifts dimension [phnx,nscreens] */
    float *yshifts)  /* array of Y fractional shifts dimension [phnx,nscreens] */
     /* ishifts[k,nscreens] and jshifts[k,nscreens] are the integer shifts for screen[k]
  xshifts[k,nscreens] and yshifts[k,nscreens] are the fractional shifts for screen[k],
     */
{
  int i,j,k,ips,jps;
  int firstel;
  float wx1,wx2,wy1,wy2;

  /* Loop on phase screens */
  for (k=0;k<nscreens;++k) {

    /* first indice for this screen */
    firstel = k*(psnx*psny);
    
    /* Loop on indices of output (integrated) phase */
    for (j=0;j<phny;++j) {
      for (i=0;i<phnx;++i) {

  ips = ishifts[i+k*phnx];
  jps = jshifts[j+k*phny];

  /* Computes the weights for the 4 surrounding pixels */
  wx1 = (1.0f - xshifts[i+k*phnx]);
  wx2 = (xshifts[i+k*phnx]);
  wy1 = (1.0f - yshifts[j+k*phny]);
  wy2 = (yshifts[j+k*phny]);


  /* Safety net: don't access elements outside of pscreens memory space */
  if ( (firstel+ips+1+(jps+1)*psnx) >= (psnx*psny*nscreens) ) {return (1);}

  /* Finaly, compute and integrate outphase */
  outphase[i+j*phnx] += ( wx1*wy1*pscreens[firstel+ips+jps*psnx]
        + wx2*wy1*pscreens[firstel+ips+1+jps*psnx]
        + wx1*wy2*pscreens[firstel+ips+(jps+1)*psnx]
        + wx2*wy2*pscreens[firstel+ips+1+(jps+1)*psnx]);

      }
    }
  }
  return(0);
}

