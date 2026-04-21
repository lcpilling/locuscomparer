#' Read association summary statistics and append normalized columns.
#' The input must contain chromosome, position, and either a -log10(p-value)
#' column or a raw p-value column. If a raw p-value column is supplied (via
#' \code{pval_col}), it is converted to -log10(p-value) automatically.
#'
#' @param in_fn (string or data.frame) Path to the input file, or a data.frame.
#' @param chromosome_col (string, optional) Name of the chromosome column. Default: 'chromosome'.
#' @param position_col (string, optional) Name of the position column. Default: 'position'.
#' @param logp_col (string, optional) Name of the -log10(p-value) column.
#'   Default: '-log10 p-value'. If this column is absent but \code{pval_col}
#'   is present, the p-value is converted automatically.
#' @param pval_col (string, optional) Name of a raw p-value column to use as a
#'   fallback when \code{logp_col} is not found. Default: 'pval'.
#' @examples
#' # Using a data.frame with a raw p-value column:
#' df = data.frame(chromosome = '1', position = 1:10, pval = 10^-(1:10))
#' d1 = read_metal(df, chromosome_col = 'chromosome', position_col = 'position', pval_col = 'pval')
#' @export
read_metal=function(in_fn, chromosome_col='chromosome', position_col='position',
                    logp_col='-log10 p-value', pval_col='pval'){
    # message('Reading ', in_fn)

    if (is.character(in_fn)){

        d = read.table(in_fn, header = TRUE, stringsAsFactors = FALSE)

    } else if (is.data.frame(in_fn)){

        d = in_fn

    } else {

        stop('The argument "in_fn" must be a string or a data.frame')

    }

    coord_cols = c(chromosome_col, position_col)
    missing_coord = setdiff(coord_cols, colnames(d))
    if (length(missing_coord) > 0) {
        stop(sprintf('Missing required columns: %s', paste(missing_coord, collapse = ', ')))
    }

    if (logp_col %in% colnames(d)) {
        logp_values = as.numeric(d[[logp_col]])
    } else if (pval_col %in% colnames(d)) {
        message(sprintf('Column "%s" not found; converting "%s" to -log10(p-value).', logp_col, pval_col))
        logp_values = -log10(as.numeric(d[[pval_col]]))
    } else {
        stop(sprintf(
            'No suitable value column found. Provide either a -log10(p-value) column ("%s") or a p-value column ("%s").',
            logp_col, pval_col
        ))
    }

    d = data.frame(
        chr = d[[chromosome_col]],
        pos = d[[position_col]],
        logp = logp_values,
        stringsAsFactors = FALSE
    )
    d$chr = as.character(d$chr)
    d$pos = as.numeric(d$pos)
    d$logp = as.numeric(d$logp)
    d$snp_id = .create_snp_id(d$chr, d$pos)
    return(d[,c('chr','pos','logp','snp_id')])
}

.create_snp_id = function(chr, pos){
    paste0(chr, ':', pos)
}

.normalize_snp_map = function(snp_map){
    if (is.null(snp_map)) return(NULL)
    required_cols = c('chromosome', 'position', 'rsid')
    missing_cols = setdiff(required_cols, colnames(snp_map))
    if (length(missing_cols) > 0) {
        stop(sprintf('The "snp" data.frame must contain columns: %s', paste(required_cols, collapse = ', ')))
    }
    res = snp_map[, required_cols]
    colnames(res) = c('chr', 'pos', 'rsid')
    res$chr = as.character(res$chr)
    res$pos = as.numeric(res$pos)
    res$snp_id = .create_snp_id(res$chr, res$pos)
    return(res)
}

.normalize_lead_ld = function(lead_ld){
    if (is.null(lead_ld)) return(NULL)
    required_cols = c('chromosome', 'position', 'r2')
    missing_cols = setdiff(required_cols, colnames(lead_ld))
    if (length(missing_cols) > 0) {
        stop(sprintf('The "lead_ld" data.frame must contain columns: %s', paste(required_cols, collapse = ', ')))
    }
    res = lead_ld[, required_cols]
    colnames(res) = c('chr', 'pos', 'r2')
    res$chr = as.character(res$chr)
    res$pos = as.numeric(res$pos)
    res$r2 = as.numeric(res$r2)
    res$snp_id = .create_snp_id(res$chr, res$pos)
    return(res)
}

#' Append two columns, chromosome (chr) and position (pos), to the input data.frame.
#'
#' @param x (data.frame) Input data.frame.
#' @param genome (string, optional) Genome assembly, either 'hg19' or 'hg38'. Default: 'hg19'.
#' @examples
#' in_fn = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn, marker_col = 'rsid', pval_col = 'pval')
#' get_position(d1, genome)
#' @export
get_position=function(x, genome = c('hg19','hg38')){

    data(config)
    on.exit(rm(config))

    conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
    on.exit(RMySQL::dbDisconnect(conn))

    stopifnot('rsid' %in% colnames(x))

    genome = match.arg(genome)

    cmd = sprintf("select rsid, chr, pos from tkg_p3v5a_%s where rsid in ('%s')",genome,paste0(x$rsid,collapse="','"))
    res = DBI::dbGetQuery(conn = conn, statement = cmd)
    y=merge(x,res,by='rsid')
    return(y)
}


#' Retrive SNP pairwise LD from database.
#' SNP pairwise lD are calculated based on 1000 Genomes Project Phase 3 version 5.
#' For storage-efficiency, the output will only include SNPs with r2 > 0.2 with the
#' input SNP.
#' @param chr (string) Chromosome name. e.g. '22'. Notice that the name should not contain 'chr'.
#' @param snp (string) SNP rsID.
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @examples
#' retrieve_LD('6', 'rs9349379', 'AFR')
#'
#' @export
retrieve_LD = function(chr,snp,population){
    data(config)
    on.exit(rm(config))

    conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
    on.exit(RMySQL::dbDisconnect(conn))

    res1 = DBI::dbGetQuery(
        conn = conn,
        statement = sprintf(
            "select SNP_A, SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_A = '%s';",
            chr,
            population,
            snp
        )
    )

    res2 = DBI::dbGetQuery(
        conn = conn,
        statement = sprintf(
            "select SNP_B as SNP_A, SNP_A as SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_B = '%s';",
            chr,
            population,
            snp
        )
    )

    res = rbind(res1,res2)
    return(res)
}

#' Get the lead SNP from the list of SNPs in input data.frame
#' The lead SNP is defined as the SNP with the lowest sum of p-values
#' from the two studies.
#' @param merged (data.frame) Input data.frame, which is a result by merging two association studies.
#' @param snp (string, optional) Lead SNP coordinate (CHR:POS). If NULL, the function will select the
#' lead SNP based on the largest sum of -log10(p-values) from the two studies.
#' @examples
#' # Select the lead SNP
#' in_fn_1 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_1, marker_col = 'rsid', pval_col = 'pval')
#' in_fn_2 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_2, marker_col = 'rsid', pval_col = 'pval')
#' merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
#' get_lead_snp(merged)
#' @export
get_lead_snp = function(merged, snp = NULL){
    lead_snp_id = snp
    if (is.null(lead_snp_id)) {
        lead_snp_id = merged[which.max(merged$logp1 + merged$logp2), 'snp_id']
    }
    else {
        if (!lead_snp_id %in% merged$snp_id) {
            stop(sprintf("%s not found in the intersection of in_fn1 and in_fn2.", lead_snp_id))
        }
    }
    return(as.character(lead_snp_id))
}

#' Assign color to each SNP according to LD.
#' @param snp_id (character vector) A vector of SNP identifiers in "CHR:POS" format on which to assign color.
#' @param snp (string) Lead SNP identifier in "CHR:POS" format. This SNP will be colored purple.
#' Other SNPs will be assigned color based on their LD with the lead SNP.
#' @param ld (data.frame) A data.frame with columns chromosome, position, and r2.
#' @examples
#' # the data.frame merged comes from the example for `get_lead_snp()`.
#' # the data.frame ld comes from the example for `retrieve_LD()`.
#' color = assign_color(snp_id = merged$snp_id, snp = '1:12345', ld)
#' @export
assign_color=function(snp_id,snp,ld=NULL){

    color = data.frame(snp_id = snp_id, stringsAsFactors = FALSE)
    color$color = 'blue4'

    ld_norm = .normalize_lead_ld(ld)
    if (!is.null(ld_norm) && nrow(ld_norm) > 0) {
        ld_norm$color = as.character(cut(ld_norm$r2, breaks=c(0,0.2,0.4,0.6,0.8,1),
                                        labels=c('blue4','skyblue','darkgreen','orange','red'),
                                        include.lowest=TRUE))
        idx = match(color$snp_id, ld_norm$snp_id)
        has_match = !is.na(idx)
        color$color[has_match] = ld_norm$color[idx[has_match]]
    }

    if (snp %in% color$snp_id){
        color[color$snp_id == snp, 'color'] = 'purple'
    }

    res = color$color
    names(res) = color$snp_id

    return(res)
}


#' Add a column of SNP labels to input data.frame
#' @param merged (data.frame) Input data.frame, which is a result by merging two
#' association studies. See the example under `get_lead_snp()` for generation of
#' such data.frame.
#' @param snp (character vector) A vector of SNP rsIDs. If only labeling one SNP,
#' this can also be a single string.
#' @examples
#' # The data.frame merged comes from the example for `get_lead_snp()`.
#' merged = add_label(merged, '1:12345')
add_label = function(merged, snp, snp_map = NULL){
    label_value = snp
    snp_map = .normalize_snp_map(snp_map)
    if (!is.null(snp_map)) {
        i = match(snp, snp_map$snp_id)
        if (!is.na(i) && !is.na(snp_map$rsid[i]) && nzchar(snp_map$rsid[i])) {
            label_value = snp_map$rsid[i]
        }
    }
    merged$label = ifelse(merged$snp_id %in% snp, label_value, '')
    return(merged)
}


#' Make a scatter plot (called the LocusCompare plot).
#' Each axis of the LocusCompare plot represent the -log10(p-value) from
#' an association study. Each point thus represent a SNP. By default, the lead SNP
#' is a purple diamond, whereas the other SNPs are colored according to
#' their LD with the lead SNP.
#' @import ggplot2
#' @import cowplot
#' @param merged (data.frame) An input data.frame which has the following
#' columns: rsid, pval1 (p-value for study 1), logp1 (p-value for study 2),
#' logp1 (log p-value for study 1), logp2 (log p-value for study 2), chr, pos.
#' See the example for `get_lead_snp()` on how to generate this data.frame.
#' @param title1 (string) The title for the x-axis.
#' @param title2 (string) The title for the y-axis.
#' @param color (data.frame) The output from `assign_color()`.
#' @param shape (data.frame) Specification of the shape of each SNP. See example blow on how to generate this data.frame.
#' @param size (data.frame) Specification of the size of each SNP. See example below on how to generate this data.frame.
#' @param legend (boolean) Whether to include the legend.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @examples
#' # The data.frame `merged` comes from the example of `add_label()`.
#' # The data.frame `color` comes from the example of `assign_color()`.
#' snp = 'rs9349379'
#' shape = ifelse(merged$rsid == snp, 23, 21)
#' names(shape) = merged$rsid
#' size = ifelse(merged$rsid == snp, 3, 2)
#' names(size) = merged$rsid
#' make_scatterplot(merged, title1 = 'GWAS', title2 = 'eQTL', color, shape, size)
#' @export
make_scatterplot = function (merged, title1, title2, color, shape, size, legend = TRUE, legend_position = c('bottomright','topright','topleft')) {

    p = ggplot(merged, aes(x = logp1, y = logp2)) +
        geom_point(aes(fill = snp_id, size = snp_id, shape = snp_id), alpha = 0.8) +
        geom_point(data = merged[merged$label != "",],
                   aes(x = logp1, y = logp2, fill = snp_id, size = snp_id, shape = snp_id)) +
        xlab(bquote(.(title1) ~ -log[10] * '(P)')) +
        ylab(bquote(.(title2) ~ -log[10] * '(P)')) +
        scale_fill_manual(values = color, guide = "none") +
        scale_shape_manual(values = shape, guide = "none") +
        scale_size_manual(values = size, guide = "none") +
        ggrepel::geom_text_repel(aes(label = label))+
        theme_classic()

    if (legend == TRUE) {
        legend_position = match.arg(legend_position)
        if (legend_position == 'bottomright'){
            legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
        } else if (legend_position == 'topright'){
            legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
        } else {
            legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
        }

        p = ggdraw(p) +
            geom_rect(data = legend_box,
                      aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                      color = "black",
                      fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
            draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
            draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
            draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
            draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
            draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 10)
    }

    return(p)
}

#' Make a locuszoom plot.
#' Details see http://locuszoom.org/.
#' @import ggplot2
#' @import cowplot
#' @param metal (data.frame) input file with two column, rsID and p-value.
#' See `read_metal()` for more details.
#' @param title (string) y-axis title.
#' @param chr (string) chromosome.
#' @param color (data.frame) The output from `assign_color()`.
#' @param shape (data.frame) Specification of the shape of each SNP. See example blow on how to generate this data.frame.
#' @param size (data.frame) Specification of the size of each SNP. See example below on how to generate this data.frame.
#' @param ylab_linebreak (boolean, optional) Whether to break the line of y-axis. If FALSE, the y-axis title and '-log10(p-value)'
#' will be on the same line. Default: FALSE.
#' @examples
#' # The data.frame `d1` comes from the example of `read_metal()`,
#' # The data.frame `color` comes from the example of `assign_color()`.
#' snp = 'rs9349379'
#' shape = ifelse(merged$rsid == snp, 23, 21)
#' names(shape) = merged$rsid
#' size = ifelse(merged$rsid == snp, 3, 2)
#' names(size) = merged$rsid
#' chr = '6'
#' make_locuszoom(d1, title = 'GWAS', chr, color, shape, size)
#' @export
make_locuszoom=function(metal,title,chr,color,shape,size,ylab_linebreak=FALSE){

    p = ggplot(metal,aes(x=pos, y=logp))+
        geom_point(aes(fill=snp_id,size=snp_id,shape=snp_id),alpha=0.8)+
        geom_point(data=metal[metal$label!='',],aes(x=pos, y=logp, fill=snp_id,size=snp_id,shape=snp_id))+
        scale_fill_manual(values=color,guide='none')+
        scale_shape_manual(values=shape,guide='none')+
        scale_size_manual(values=size,guide='none')+
        scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)})+
        ggrepel::geom_text_repel(aes(label=label))+
        xlab(paste0('chr',chr,' (Mb)'))+
        ylab(bquote(.(title)~-log[10]*'(P)'))+
        theme_classic()+
        theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"))

    if (ylab_linebreak==TRUE){
        p = p + ylab(bquote(atop(.(title),-log[10]*'(P)')))
    }
    return(p)
}


#' Generated a combined plot with two locuszoom plots and a locuscompare
#' plot. Each locuszoom plot represent an association study.
#' @param merged (data.frame) An input data.frame which has the following
#' columns: rsid, pval1 (p-value for study 1), logp1 (p-value for study 2),
#' logp1 (log p-value for study 1), logp2 (log p-value for study 2), chr, pos.
#' See the example for `get_lead_snp()` on how to generate this data.frame.
#' @param title1 (string) The title for the x-axis.
#' @param title2 (string) The title for the y-axis.
#' @param ld (data.frame) The output from `retrieve_LD()`.
#' @param chr (string) Chromosome name. e.g. '22'. Notice that the name should not contain 'chr'.
#' @param snp (string, optional) SNP rsID. If NULL, the function will select the lead SNP. Default: NULL.
#' @param combine (boolean, optional) Should the three plots be combined into one plot? If FALSE, a list of
#' three plots will be returned. Default: TRUE.
#' @param legend (boolean, optional) Should the legend be shown? Default: TRUE.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @param lz_ylab_linebreak (boolean, optional) Whether to break the line of y-axis of the locuszoom plot.
#' If FALSE, the y-axis title and '-log10(p-value)'. will be on the same line. Default: FALSE.
#' @examples
#' # The data.frame `merged` comes from the example of `add_label()`.
#' # the data.frame `ld` comes from the example for `retrieve_LD()`.
#' make_combined_plot(merged, 'GWAS', 'eQTL', ld, chr)
#' @export
make_combined_plot = function (merged, title1, title2, ld, chr, snp = NULL, snp_map = NULL,
                               combine = TRUE, legend = TRUE,
                               legend_position = c('bottomright','topright','topleft'),
                               lz_ylab_linebreak=FALSE) {

    snp = get_lead_snp(merged, snp)
    # print(sprintf("INFO - %s", snp))

    color = assign_color(merged$snp_id, snp, ld)

    shape = ifelse(merged$snp_id == snp, 23, 21)
    names(shape) = merged$snp_id

    size = ifelse(merged$snp_id == snp, 3, 2)
    names(size) = merged$snp_id

    merged = add_label(merged, snp, snp_map)

    p1 = make_scatterplot(merged, title1, title2, color,
                          shape, size, legend, legend_position)

    metal1 = merged[,c('snp_id', 'logp1', 'chr', 'pos', 'label')]
    colnames(metal1)[which(colnames(metal1) == 'logp1')] = 'logp'
    p2 = make_locuszoom(metal1, title1, chr, color, shape, size, lz_ylab_linebreak)

    metal2 = merged[,c('snp_id', 'logp2', 'chr', 'pos', 'label')]
    colnames(metal2)[which(colnames(metal2) == 'logp2')] = 'logp'
    p3 = make_locuszoom(metal2, title2, chr, color, shape, size, lz_ylab_linebreak)

    if (combine) {
        p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
        p4 = cowplot::plot_grid(p2, p3, align = "v", nrow = 2, rel_heights=c(0.8,1))
        p5 = cowplot::plot_grid(p1, p4)
        return(p5)
    }
    else {
        return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
    }
}

#' Make a locuscompare plot.
#' @param in_fn1 (string or data.frame) Path to the input file for study 1, or a data.frame.
#' @param in_fn2 (string or data.frame) Path to the input file for study 2, or a data.frame.
#' @param chromosome_col1 (string, optional) Name of the chromosome column in dataset 1. Default: 'chromosome'.
#' @param position_col1 (string, optional) Name of the position column in dataset 1. Default: 'position'.
#' @param logp_col1 (string, optional) Name of the -log10(p-value) column in dataset 1. Default: '-log10 p-value'.
#'   If absent, \code{pval_col1} is used and converted automatically.
#' @param pval_col1 (string, optional) Fallback raw p-value column for dataset 1. Default: 'pval'.
#' @param title1 (string) The title for the x-axis.
#' @param chromosome_col2 (string, optional) Name of the chromosome column in dataset 2. Default: 'chromosome'.
#' @param position_col2 (string, optional) Name of the position column in dataset 2. Default: 'position'.
#' @param logp_col2 (string, optional) Name of the -log10(p-value) column in dataset 2. Default: '-log10 p-value'.
#'   If absent, \code{pval_col2} is used and converted automatically.
#' @param pval_col2 (string, optional) Fallback raw p-value column for dataset 2. Default: 'pval'.
#' @param title2 (string) The title for the y-axis.
#' @param snp (string or data.frame, optional) Either a lead SNP identifier ("CHR:POS") or a
#' data.frame with columns chromosome, position, and rsid for optional RSID labeling. Default: NULL.
#' @param lead_ld (data.frame, optional) A data.frame with columns chromosome, position, and r2 for
#' LD with the lead SNP. If omitted, a fallback database retrieval is attempted and silently ignored if unavailable.
#' @param population (string, optional) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'. Default: 'EUR'.
#' @param min_match (integer, optional) Minimum number of overlapping variants required between datasets.
#' If fewer overlaps are found, an error is raised. Default: 10.
#' @param combine (boolean, optional) Should the three plots be combined into one plot? If FALSE, a list of
#' three plots will be returned. Default: TRUE.
#' @param legend (boolean, optional) Should the legend be shown? Default: TRUE.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @param lz_ylab_linebreak (boolean, optional) Whether to break the line of y-axis of the locuszoom plot.
#' @examples
#' # Using data.frames with a raw p-value column:
#' d1 = data.frame(chromosome = '1', position = 1:100, pval = 10^-runif(100, 1, 10))
#' d2 = data.frame(chromosome = '1', position = 1:100, pval = 10^-runif(100, 1, 10))
#' locuscompare(in_fn1 = d1, in_fn2 = d2, min_match = 10)
#' @export
locuscompare = function(in_fn1, in_fn2, chromosome_col1 = "chromosome", position_col1 = "position",
                 logp_col1 = "-log10 p-value", pval_col1 = "pval", title1 = "eQTL",
                 chromosome_col2 = "chromosome", position_col2 = "position",
                 logp_col2 = "-log10 p-value", pval_col2 = "pval", title2 = "GWAS",
                 snp = NULL, lead_ld = NULL, population = "EUR", min_match = 10, combine = TRUE, legend = TRUE,
                 legend_position = c('bottomright','topright','topleft'),
                 lz_ylab_linebreak = FALSE) {
    d1 = read_metal(in_fn1, chromosome_col1, position_col1, logp_col1, pval_col1)
    d2 = read_metal(in_fn2, chromosome_col2, position_col2, logp_col2, pval_col2)

    merged = merge(d1, d2, by = c("chr", "pos", "snp_id"), suffixes = c("1", "2"), all = FALSE)
    if (nrow(merged) < min_match) {
        stop(sprintf('Only %d overlapping variants were found between in_fn1 and in_fn2 (minimum required: %d); possible causes include genome build mismatch or dataset filtering differences.', nrow(merged), min_match))
    }

    chr = unique(merged$chr)
    if (length(chr) != 1) stop('There must be one and only one chromosome.')

    snp_map = NULL
    lead_snp = NULL
    if (!is.null(snp)) {
        if (is.data.frame(snp)) {
            snp_map = .normalize_snp_map(snp)
        } else if (is.character(snp) && length(snp) == 1) {
            lead_snp = snp
        } else {
            stop('The "snp" argument must be NULL, a lead SNP string (CHR:POS), or a data.frame with chromosome, position, rsid.')
        }
    }

    lead_snp = get_lead_snp(merged, lead_snp)

    ld = .normalize_lead_ld(lead_ld)
    if (is.null(ld) && !is.null(snp_map)) {
        lead_idx = match(lead_snp, snp_map$snp_id)
        if (!is.na(lead_idx)) {
            lead_rsid = snp_map$rsid[lead_idx]
        } else {
            lead_rsid = NA
        }
        if (!is.na(lead_rsid) && nzchar(lead_rsid)) {
            ld_backup = tryCatch(
                retrieve_LD(chr, lead_rsid, population),
                error = function(e) NULL
            )
            if (!is.null(ld_backup) && nrow(ld_backup) > 0) {
                snp_lookup = snp_map[, c('snp_id', 'rsid', 'chr', 'pos')]
                ld_join = merge(ld_backup[, c('SNP_B', 'R2')], snp_lookup, by.x = 'SNP_B', by.y = 'rsid', all.x = FALSE)
                if (nrow(ld_join) > 0) {
                    ld = .normalize_lead_ld(data.frame(chromosome = ld_join$chr, position = ld_join$pos, r2 = ld_join$R2, stringsAsFactors = FALSE))
                }
            }
        }
    }

    p = make_combined_plot(merged, title1, title2, ld, chr, lead_snp, snp_map, combine,
                           legend, legend_position, lz_ylab_linebreak)
    return(p)
}
