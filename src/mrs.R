mrs <- function(nonreference, reference){
	nonreference.order = ceiling(rank(nonreference))
	reference[nonreference.order]
}
