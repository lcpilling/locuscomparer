context('Test assign color')

positions = 0:11
rsid = paste0('1:', positions)
ld = data.frame(
    chromosome = '1',
    position = 0:10,
    r2 = seq(1, 0, -0.1)
)

res = assign_color(rsid, '1:1', ld)

test_that('assign_color',{
    expect_equal(unname(res['1:1']),'purple')
    expect_equal(unname(res['1:11']),'blue4')
    expect_equal(unname(res['1:0']),'red')
})
