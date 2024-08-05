#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})
suppressMessages(library(knitr))
suppressMessages(library(rmarkdown))

knitr::opts_chunk$set(
    warning=TRUE,
    message=TRUE,
    cache=FALSE,
    comment="",
    fig.show="asis",
    dev="png"
)

work_dir <- getwd()
temp_dir <- tempdir()
script_location <- commandArgs(trailingOnly=TRUE)[1]
knit_md_location <- file.path(
    temp_dir,
    gsub(".r$|.R$", ".knit.md", basename(script_location))
)
print(paste("Running", script_location, "through the HTML report wrapper"))
print(paste("Working directory:", work_dir))
print(paste("Temporary directory for generating HTML report:", temp_dir))

tryCatch(
    expr = {
        rmarkdown::render(
            script_location,
            output_format="html_document",
            output_file="sc_report.html",
            output_dir=work_dir,
            intermediates_dir=temp_dir,
            knit_root_dir=work_dir,
            runtime="static",
            quiet=TRUE,
            clean=FALSE
        )
    },
    error = function(e){
        print(paste("Failed to run", script_location, "with error -", e))
        rmarkdown::render(
            knit_md_location,
            output_format="html_document",
            output_file="sc_report.html",
            output_dir=work_dir,
            output_options=list(
                pandoc_args=c(
                    "--metadata", paste0("title=", basename(script_location))
                )
            ),
            intermediates_dir=temp_dir,
            knit_root_dir=work_dir,
            runtime="static",
            quiet=TRUE,
            clean=FALSE
        )
        quit(save="no", status=1, runLast=FALSE)
    }
)

