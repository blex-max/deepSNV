#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "Rinternals.h"
#include "deepsnv-prototypes.h"

static R_NativePrimitiveArgType C_bam2R_type[14] = {
  STRSXP, STRSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static const R_CMethodDef CEntries[] = {
    {"C_bam2R",       (DL_FUNC) &bam2R,       14, C_bam2R_type},
    {"C_dbetabinom",  (DL_FUNC) &dbetabinom,  11},
    {"C_pbetabinom",  (DL_FUNC) &pbetabinom,  11},
    {NULL, NULL, 0}
};

void R_init_deepSNV(DllInfo *info) {
    R_registerRoutines(info, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
