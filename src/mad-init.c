#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bDeviation(char **, double *, int *, double *);

static const R_CMethodDef CEntries[] = {
  {"bDeviation",               (DL_FUNC) &bDeviation,               4},
  {NULL, NULL, 0}
};

void R_init_mad(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}