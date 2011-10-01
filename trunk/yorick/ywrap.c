/* codger-generated yorick package wrapper file */
#include "play.h"
#include "ydata.h"

/*----------------begin pray_utils.i */
extern BuiltIn Y__get2dPhase;

extern int _get2dPhase(void *, int , int , int , void *, int , 
  int , void *, void *, void *, void *);
void
Y__get2dPhase(int n)
{
  if (n!=11) YError("_get2dPhase takes exactly 11 arguments");
  PushIntValue(_get2dPhase(yarg_sp(10), yarg_si(9), yarg_si(8), 
    yarg_si(7), yarg_sp(6), yarg_si(5), yarg_si(4), yarg_sp(3), 
    yarg_sp(2), yarg_sp(1), yarg_sp(0)));
}

extern BuiltIn Y__dist;

extern int _dist(void *, long , long , float , float );
void
Y__dist(int n)
{
  if (n!=5) YError("_dist takes exactly 5 arguments");
  PushIntValue(_dist(yarg_sp(4), yarg_sl(3), yarg_sl(2), 
    yarg_sf(1), yarg_sf(0)));
}


/*----------------list include files */

static char *y0_includes[] = {
  "pray_utils.i",
  0,
  0
};

/*----------------collect pointers and names */

static BuiltIn *y0_routines[] = {
  &Y__get2dPhase,
  &Y__dist,
  0
};

static void *y0_values[] = {
  0
};

static char *y0_names[] = {
  "_get2dPhase",
  "_dist",
  0
};

/*----------------define package initialization function */

PLUG_EXPORT char *yk_pray_utils(char ***,
                         BuiltIn ***, void ***, char ***);
static char *y0_pkgname = "pray_utils";

char *
yk_pray_utils(char ***ifiles,
       BuiltIn ***code, void ***data, char ***varname)
{
  *ifiles = y0_includes;
  *code = y0_routines;
  *data = y0_values;
  *varname = y0_names;
  return y0_pkgname;
}
