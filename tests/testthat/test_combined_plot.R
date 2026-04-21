context('Test combined plot')

d1 = data.frame(chromosome = '1', position = 1:12, check.names = FALSE)
d1[['-log10 p-value']] = seq(1, 12)
d2 = data.frame(chromosome = '1', position = 1:12, check.names = FALSE)
d2[['-log10 p-value']] = seq(12, 1)
rd1 = read_metal(d1)
rd2 = read_metal(d2)
merged = merge(rd1, rd2, by = c("chr", "pos", "snp_id"), suffixes = c("1", "2"), all = FALSE)
chr = unique(merged$chr)
ld = data.frame(chromosome = '1', position = 1:12, r2 = seq(0, 1, length.out = 12))

test_that('make_combined_plot returns a combined plot object',{
    p = make_combined_plot(merged, 'GWAS', 'eQTL', ld, chr, snp = NULL,
                           combine = TRUE, legend = TRUE,
                           legend_position = 'bottomright', lz_ylab_linebreak=FALSE)
    expect_true(inherits(p, 'gg') || inherits(p, 'gtable') || inherits(p, 'ggplot'))
})
