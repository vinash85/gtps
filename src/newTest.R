calculate.r21 <- function(n, d){
	N = sum(n) 
	D = sum(d)
	t1 = (n[1] + n[2]) * d[2]/n[2] + (n[4] + n[3]) * d[4]/n[4] - D
	(D - (n[4] + n[2]) * t1/ (N-D))/(n[1] + n[3]) 
}

interactiontest <- function(n, events, quadrant) {
# events: 0 - death; 1 - censored
r21 = calculate.r21(n, d)
num =  ( d[3] - r21 * n[3] ) 
den = n[3] * r21 * (1 - r21)
for (tt in seq(length(events)) {
	event = events[tt]
	if(event){ # censored 
	n[quadrant[tt]]  = n[quadrant[tt]] - 1   
	} else{ # death
		d[quadrant[tt]] = d[quadrant[tt]] + 1
	}
	r21 = calculate.r21(n, d)
	num =  num + ( d[3] - r21 * n[3] ) 
	den = den + n[3] * r21 * (1 - r21)
	}
	zscore = num/sqrt(den)
	p.value = dnorm(zscore)
	list(zscore, p.value) 
}