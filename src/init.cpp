#include "DropletUtils.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(downsample_matrix, 3),
    REGISTER(downsample_runs, 4),
    REGISTER(montecarlo_pval, 5),
    REGISTER(compute_multinom, 2),
    REGISTER(find_swapped, 3),
    REGISTER(get_cell_barcodes, 3),
    REGISTER(load_tenx_to_hdf5, 6),
    {NULL, NULL, 0}
};

void attribute_visible R_init_DropletUtils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}
