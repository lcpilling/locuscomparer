context('test lead SNP selection')

merged = data.frame(
    chr = rep('1', 3),
    pos = c(100, 200, 300),
    snp_id = c('1:100', '1:200', '1:300'),
    logp1 = c(2, 6, 3),
    logp2 = c(1, 4, 10),
    stringsAsFactors = FALSE
)

test_that('get_lead_snp uses max sum of logp values',{
    expect_equal(get_lead_snp(merged), '1:300')
})
