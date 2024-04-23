import("future", attach=FALSE)
import("data.table", attach=FALSE)
import("BiocParallel", attach=FALSE)

export(
    "parallel"
)


parallel <- function (args) {
    base::print(
        base::paste(
            "Setting parallelization to", args$cpus, "cores, and", args$memory,
            "GB of memory allowed to be shared between the processes"
        )
    )
    invisible(utils::capture.output(future::plan("multiprocess", workers=args$cpus)))
    invisible(utils::capture.output(future::plan()))
    invisible(utils::capture.output(data.table::setDTthreads(args$cpus)))
    base::options(future.globals.maxSize = args$memory * 1024^3)                        # convert to bytes
    BiocParallel::register(BiocParallel::MulticoreParam(args$cpus, RNGseed=args$seed))  # for DESeq2, RNGseed is hardcoded to make results reproducible
    base::set.seed(args$seed)
}