context('Test read metal')

in_fn1 = data.frame(
    chromosome = c('1', '1'),
    position = c(100, 200),
    check.names = FALSE
)
in_fn1[['-log10 p-value']] = c(2, 3)
d1 = read_metal(in_fn1, chromosome_col = 'chromosome', position_col = 'position', logp_col = '-log10 p-value')

test_that('headers are correct',{
    expect_equal(colnames(d1), c('chr','pos','logp','snp_id'))
    expect_equal(d1$snp_id, c('1:100', '1:200'))
})
