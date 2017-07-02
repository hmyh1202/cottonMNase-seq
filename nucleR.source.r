df2gff <- function (df, ...)
{   # Convert a data.frame into a form that is easily saved in disk as a gff.
    # It tries to match gff fields to the row names of the data.frame.
    # gff fields not found in the data.frame will be set to `.` and fields in
    # the data.frame not found in the gff fields will be stored accordingly in
    # the `attribute` field, separated by `;`.
    # It also accepts optional keyword arguments to specify gff fields that
    # are not present in the data.frame.
    kwargs <- list(...)
    fields <- c("seqname", "source", "feature", "start", "end", "score",
    "strand", "frame")
    out.df <- data.frame(matrix(nrow=nrow(df),
    ncol=length(fields) + 1,
    dimnames=list(c(),
    c(fields, "attribute"))))
    for (f in fields) {
        if (f %in% colnames(df)) {
            out.df[[f]] <- as.vector(df[[f]])
        } else if (f %in% names(kwargs)) {
            if (length(kwargs[[f]]) == 1) {
                out.df[[f]] <- rep(kwargs[[f]], nrow(out.df))
            } else {
                out.df[[f]] <- kwargs[[f]]
            }
        } else {
            out.df[[f]] <- rep(".", nrow(out.df))
        }
    }
    nonfield.columns <- colnames(df)[!colnames(df) %in% fields]
    attrVal <- function(i, df) sprintf("%s=%s", i, df[[i]])
    attrs <- do.call(paste,
    c(lapply(nonfield.columns,
    attrVal,
    df),
    sep=";"))
    if (length(attrs)) {
        out.df[["attribute"]] <- attrs
    } else {
        out.df[["attribute"]] <- rep(".", nrow(out.df))
    }
    out.df
}

writeGff <- function (df, outpath)
{
    # Use this to write the output of df2gff to disk.
    write.table(df,
    file=outpath,
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE)
}
