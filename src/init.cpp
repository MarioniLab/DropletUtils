#include "DropletUtils.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(downsample_matrix, 2),
    REGISTER(calculate_random_dev, 2),
    REGISTER(calculate_pval, 5),
    {NULL, NULL, 0}
};

void attribute_visible R_init_DropletUtils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}
