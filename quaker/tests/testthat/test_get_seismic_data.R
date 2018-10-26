library(quaker)

context('Checks that get_seismic_data works as expected')

test_that(desc = 'get_seismic_data produces geojson ', code = {
  expect_output(object = str(get_seismic_data(timeFrame = 'PAST_DAY', minMagnitude = '1')), 'List of 4')
  expect_is(object = get_seismic_data(timeFrame = 'PAST_WEEK', minMagnitude = 'all'), class = 'seismic_geojson')
  expect_error(get_seismic_data(timeFrame = 'PAST_CENTUARY', minMagnitude = 'all'))
})

