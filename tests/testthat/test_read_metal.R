context('Test read metal')

# Test using explicit -log10 p-value column
in_fn1 = data.frame(
    chromosome = c('1', '1'),
    position = c(100, 200),
    check.names = FALSE
)
in_fn1[['-log10 p-value']] = c(2, 3)
d1 = read_metal(in_fn1, chromosome_col = 'chromosome', position_col = 'position', logp_col = '-log10 p-value')

test_that('headers are correct with logp column',{
    expect_equal(colnames(d1), c('chr','pos','logp','snp_id'))
    expect_equal(d1$snp_id, c('1:100', '1:200'))
    expect_equal(d1$logp, c(2, 3))
})

# Test pval fallback: if only 'pval' column is present, it should be auto-converted
in_fn2 = data.frame(
    chromosome = c('1', '1'),
    position = c(100, 200),
    pval = c(0.01, 0.001)
)
d2 = suppressMessages(read_metal(in_fn2, chromosome_col = 'chromosome', position_col = 'position', pval_col = 'pval'))

test_that('pval is correctly converted to -log10(pval)',{
    expect_equal(colnames(d2), c('chr','pos','logp','snp_id'))
    expect_equal(d2$logp, -log10(c(0.01, 0.001)))
})

test_that('error raised when neither logp nor pval column found',{
    in_fn3 = data.frame(chromosome = '1', position = 1, other_col = 0.05)
    expect_error(read_metal(in_fn3))
})
