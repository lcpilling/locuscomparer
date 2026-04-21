context('Test locuscompare')

d1 = data.frame(chromosome = '1', position = 1:12, check.names = FALSE)
d1[['-log10 p-value']] = seq(1, 12)
d2 = data.frame(chromosome = '1', position = 1:12, check.names = FALSE)
d2[['-log10 p-value']] = seq(12, 1)
lead_ld = data.frame(chromosome = '1', position = 1:12, r2 = seq(0, 1, length.out = 12))
snp_map = data.frame(chromosome = '1', position = 12, rsid = 'rsLead', stringsAsFactors = FALSE)

test_that('locuscompare works with chr/pos/logp input and optional lead_ld',{
    p = locuscompare(in_fn1 = d1, in_fn2 = d2, lead_ld = lead_ld, snp = snp_map, min_match = 10)
    expect_true(inherits(p, 'gg') || inherits(p, 'gtable') || inherits(p, 'ggplot'))
})

test_that('locuscompare works with raw pval columns (auto-converts to -log10)',{
    d1_pval = data.frame(chromosome = '1', position = 1:12, pval = 10^-seq(1, 12))
    d2_pval = data.frame(chromosome = '1', position = 1:12, pval = 10^-seq(12, 1))
    p = suppressMessages(locuscompare(in_fn1 = d1_pval, in_fn2 = d2_pval, min_match = 10))
    expect_true(inherits(p, 'gg') || inherits(p, 'gtable') || inherits(p, 'ggplot'))
})

test_that('locuscompare errors on too few overlaps',{
    d2_small = d2[d2$position <= 2, ]
    expect_error(
        locuscompare(in_fn1 = d1, in_fn2 = d2_small, min_match = 10),
        regexp = 'overlapping variants'
    )
})
