safeBPPARAM <- if (.Platform$OS.type=="windows") BiocParallel::SerialParam() else BiocParallel::MulticoreParam(3)
