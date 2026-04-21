context('Test locuscatter')

d1 = data.frame(chromosome = '1', position = 1:12, check.names = FALSE)
d1[['-log10 p-value']] = seq(1, 12)
d2 = data.frame(chromosome = '1', position = 1:12, check.names = FALSE)
d2[['-log10 p-value']] = seq(12, 1)

rd1 = read_metal(d1)
rd2 = read_metal(d2)
merged = merge(rd1, rd2, by = c("chr", "pos", "snp_id"), suffixes = c("1", "2"), all = FALSE)
snp = get_lead_snp(merged)
lead_ld = data.frame(chromosome = '1', position = 1:12, r2 = seq(0, 1, length.out = 12))
color = assign_color(merged$snp_id, snp, lead_ld)

shape = ifelse(merged$snp_id == snp, 23, 21)
names(shape) = merged$snp_id

size = ifelse(merged$snp_id == snp, 3, 2)
names(size) = merged$snp_id

merged = add_label(merged, snp)

test_that('make_scatterplot returns a ggplot object',{
    p = make_scatterplot(merged, title1 = 'GWAS', title2 = 'eQTL', color, shape,
                         size, legend = TRUE, legend_position = 'bottomright')
    expect_true(inherits(p, 'gg'))
})
