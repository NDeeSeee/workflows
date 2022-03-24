import("dplyr", attach=FALSE)
import("purrr", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("tibble", attach=FALSE)
import("SeuratDisk", attach=FALSE)
import("rtracklayer", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)

export(
    "get_file_type",
    "load_barcodes_data",
    "export_data",
    "export_rds",
    "load_cell_identity_data",
    "extend_metadata",
    "load_grouping_data",
    "load_blacklisted_data",
    "load_10x_multiome_data"
)


get_file_type <- function (filename){
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}

load_barcodes_data <- function(location, seurat_data){
    default_barcodes_data <- SeuratObject::Cells(seurat_data)                # to include all available cells
    if (!is.null(location)){
        barcodes_data <- utils::read.table(
            location,
            sep=get_file_type(location),
            header=FALSE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )[,1]                                                                 # to get it as vector
        base::print(
            base::paste(
                "Barcodes data is successfully loaded and from", location,
                "and will be used to prefilter feature-barcode matrices by",
                "cells of interest")
        )
        return (barcodes_data)
    }
    base::print("Barcodes data is not provided. Using all cells")
    return (default_barcodes_data)
}

export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE){
    base::tryCatch(
        expr = {
            utils::write.table(
                data,
                file=location,
                sep=get_file_type(location),
                row.names=row_names,
                col.names=col_names,
                quote=quote
            )
            base::print(base::paste("Export data to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export data to", location, sep=" "))
        }
    )
}

export_rds <- function(data, location){
    base::tryCatch(
        expr = {
            base::saveRDS(data, location)
            base::print(base::paste("Export data as RDS to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export data as RDS to", location, sep=" "))
        }
    )
}

export_h5seurat <- function(data, location, overwrite=TRUE){
    base::tryCatch(
        expr = {
            SeuratDisk::SaveH5Seurat(data, location, overwrite=overwrite)
            base::print(base::paste("Export data as h5seurat to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export data as h5seurat to", location, sep=" "))
        }
    )
}

load_cell_identity_data <- function(location) {
    cell_identity_data <- utils::read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    # prepend with LETTERS, otherwise the order on the plot will be arbitrary sorted
    cell_identity_data <- cell_identity_data %>%
                          dplyr::mutate("library_id"=base::paste(LETTERS[1:base::nrow(cell_identity_data)], .$library_id))
    return (cell_identity_data)
}

extend_metadata <- function(seurat_data, location){
    base::tryCatch(
        expr = {
            extra_metadata <- utils::read.table(
                location,
                sep=get_file_type(location),
                header=TRUE,
                check.names=FALSE,
                stringsAsFactors=FALSE
            )
            base::print(base::paste("Extra metadata is successfully loaded from ", location))
            refactored_metadata <- base::data.frame(SeuratObject::Cells(seurat_data)) %>%      # create a dataframe with only one column
                                   dplyr::rename("barcode"=1) %>%                        # rename that column to barcode
                                   dplyr::left_join(extra_metadata, by="barcode") %>%    # intersect with loaded extra metadata by "barcode"
                                   tibble::remove_rownames() %>%
                                   tibble::column_to_rownames("barcode") %>%
                                   replace(is.na(.), "Unknown") %>%             # in case an extra metadata had less barcodes than we had in our Seurat object
                                   dplyr::rename_with(~base::paste0("custom_", .x))          # add prefix to all extra metadata columns
            seurat_data <- SeuratObject::AddMetaData(
                seurat_data,
                refactored_metadata[Cells(seurat_data), , drop=FALSE]           # to guarantee the proper cells order
            )
        },
        error = function(e){
            base::print(base::paste("Failed to add extra cells metadata due to", e))
        },
        finally = {
            return (seurat_data)
        }
    )
}

load_grouping_data <- function(location, cell_identity_data) {
    default_grouping_data <- base::data.frame(
        library_id=cell_identity_data$library_id,
        condition=cell_identity_data$library_id,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    if (!is.null(location)){
        grouping_data <- utils::read.table(
            location,
            sep=get_file_type(location),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        # prepend with LETTERS to correspond to the library_id from the cell_identity_data
        grouping_data <- grouping_data %>%
                          dplyr::mutate("library_id"=base::paste(LETTERS[1:base::nrow(grouping_data)], .$library_id))
        if ( (base::nrow(grouping_data) == base::nrow(cell_identity_data)) && all(base::is.element(cell_identity_data$library_id, grouping_data$library_id)) ){
            base::print(base::paste("Grouping data is successfully loaded from ", location))
            return (grouping_data)
        } else {
            base::print(base::paste("Applying defaults - grouping data loaded from", location, "is malformed"))
            return (default_grouping_data)
        }
    }
    base::print("Grouping data is not provided. Applying defaults")
    return (default_grouping_data)
}

load_blacklisted_data <- function(location) {
    default_blacklisted_data <- NULL
    if (!is.null(location)){
        blacklisted_data <- rtracklayer::import(location, format="BED")
        base::print(base::paste("Blacklisted regions data is successfully loaded from ", location))
        return (blacklisted_data)
    }
    base::print("Blacklisted regions data is not provided. Cells won't be filtered by --maxblacklisted")
    return (default_blacklisted_data)
}

load_10x_multiome_data <- function(args, cell_identity_data, grouping_data) {
    base::suppressMessages(raw_data <- Seurat::Read10X(data.dir=args$mex))
    seurat_data <- SeuratObject::CreateSeuratObject(
        counts=raw_data$`Gene Expression`,
        min.cells=args$gexmincells,
        names.delim="-",
        names.field=2
    )
    annotation <- rtracklayer::import(args$annotations, format="GFF")
    seurat_data[["ATAC"]] <- Signac::CreateChromatinAssay(
        counts=raw_data$Peaks,
        sep=c(":", "-"),
        fragments=args$fragments,
        min.cells=args$atacmincells,
        annotation=annotation
    )
    base::rm(raw_data, annotation)                                             # removing unused data
    base::print("Assigning new dataset identities")
    SeuratObject::Idents(seurat_data) <- "orig.ident"                                  # safety measure to make sure we get correct Idents
    idents <- as.numeric(as.character(SeuratObject::Idents(seurat_data)))              # need to properly convert factor to numeric vector
    new_ident <- cell_identity_data$library_id[idents]
    if (sum(is.na(new_ident)) > 0){
        base::print("Identity file includes less than expected number of rows. Exiting.")
        base::quit(save="no", status=1, runLast=FALSE)
    }
    seurat_data[["new.ident"]] <- new_ident
    seurat_data[["condition"]] <- grouping_data$condition[base::match(seurat_data$new.ident, grouping_data$library_id)]
    SeuratObject::Idents(seurat_data) <- "new.ident"
    if ( base::nrow(cell_identity_data) > length(base::unique(base::as.vector(as.character(SeuratObject::Idents(seurat_data))))) ){
        base::print("Identity file includes more than expected number of rows. Exiting.")
        base::quit(save="no", status=1, runLast=FALSE)
    }
    base::gc()
    return (seurat_data)
}