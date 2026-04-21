context('Test locuszoom')

title = 'GWAS'

d1 = data.frame(chromosome = '1', position = 1:12, check.names = FALSE)
d1[['-log10 p-value']] = seq(1, 12)
d2 = data.frame(chromosome = '1', position = 1:12, check.names = FALSE)
d2[['-log10 p-value']] = seq(12, 1)
rd1 = read_metal(d1)
rd2 = read_metal(d2)
merged = merge(rd1, rd2, by = c("chr", "pos", "snp_id"), suffixes = c("1", "2"), all = FALSE)
chr = unique(merged$chr)

snp = get_lead_snp(merged)
lead_ld = data.frame(chromosome = '1', position = 1:12, r2 = seq(0, 1, length.out = 12))
color = assign_color(merged$snp_id, snp, lead_ld)

shape = ifelse(merged$snp_id == snp, 23, 21)
names(shape) = merged$snp_id

size = ifelse(merged$snp_id == snp, 3, 2)
names(size) = merged$snp_id

merged = add_label(merged, snp)
metal = merged[, c('snp_id', 'logp1', 'chr', 'pos', 'label')]
colnames(metal)[which(colnames(metal) == 'logp1')] = 'logp'

test_that('make_locuszoom returns a ggplot object',{
    p = make_locuszoom(metal,title,chr,color,shape,size,ylab_linebreak=FALSE)
    expect_true(inherits(p, 'gg'))
})
