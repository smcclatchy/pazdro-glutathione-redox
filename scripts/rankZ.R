# There is a mistake in the rankZ function that I have been using (and sharing) for years.
# I feel stupid for not noticing before but the handling of missing values was not correct
#
# The correct rankZ function should be

###
# rankZ function definition
rankZ <- function(x) {
  x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))
}

# The older one was
# ###
# # rankZ function definition
# rankZ <- function(x) {
#   x <- rank(x, na.last = "keep") / (length(x) + 1)
#   return(qnorm(x))
# }
#
# The denominator in the step to scale ranks between zero and one was incorrect.
#
# If you have written your own – you are probably fine.  Otherwise please adopt this new version.
# It only matters for traits with missing data…so it is up to you to determine the impact.
