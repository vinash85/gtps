calculate.r01 <- function(n, d){
	N = sum(n) 
	D = sum(d)
	t1 = (n[1] + n[2]) * d[2]/n[2] + (n[4] + n[3]) * d[4]/n[4] - D
	(D - (n[4] + n[2]) * t1/ (N-D))/(n[1] + n[3])
	
	  # r21 = ((n[1] + n[2]) *(D * (d[4]/n[4]-1)-(n[2] + n[4]) * d[4]/n[4] * d[2]/n[2]+(n[2] + n[4]) * d[2]/n[2]- (n[1] + n[3]) * d[4]/n[4]+ (n[1] + n[3]) * d[2]/n[2])+(d[4]/n[4]-1) *((n[3] + n[4]) * (D-(n[2] + n[4]) * d[4]/n[4])+D * (-D+(n[2] + n[4])+(n[1] + n[3]))))/((n[1] + n[3]) *((n[1] + n[2]) * (d[2]/n[2]-1)+(n[3] + n[4]) * (d[4]/n[4]-1))) 
	  
	  rb1 = ((n[3] + n[4]) *d[4]/n[4]-(n[3] + n[4]) * d[2]/n[2]+D * d[2]/n[2]-D)/((n[1] + n[2]) * d[2]/n[2]-(n[1] + n[2])+(n[3] + n[4]) *d[4]/n[4]-(n[3] + n[4])) 
	  rb2 = (-(n[1] + n[2]) * d[4]/n[4]+(n[1] + n[2]) * d[2]/n[2]+D *d[4]/n[4]-D)/((n[1] + n[2]) *d[2]/n[2]-(n[1] + n[2])+(n[3] + n[4]) * d[4]/n[4]-(n[3] + n[4])) 
	  ra1 = (-(n[1] + n[2]) * D+(n[1] + n[2]) * (n[2] + n[4]) * d[2]/n[2]-(n[3] + n[4]) * D+(n[3] + n[4]) * (n[2] + n[4]) * d[4]/n[4]+D^2-D * (n[2] + n[4]))/((n[1] + n[3]) * (-(n[1] + n[2])-(n[3] + n[4])+D))
	  r21 = ra1 + rb2 -ra1 * rb2
	  r11 = ra1 + rb1 -ra1 * rb1
	  r21 =ifelse(r21 <1, r21, 1)
	  r11 =ifelse(r11 <1, r11, 1)
	  if((r11 ==1) |(r21==1) ) browser()
	  return(c(r11, r21))
}



interactionTest.zscore <- function(n, events, quadrant) {
# events: 0 - death; 1 - censored
d= rep(0,4)
r21 = calculate.r01(n, d)[2]
num =  ( d[3] - r21 * n[3] ) 
den = n[3] * r21 * (1 - r21)
cur = rep(0, (length(events) -1))
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC403858/pdf/bmj32801073.pdf
for (tt in seq(length(events) -1)) {
	event = events[tt]
	if(event){ # censored 
	n[quadrant[tt]]  = n[quadrant[tt]] - 1   
	} else{ # death
		d[quadrant[tt]] = d[quadrant[tt]] + 1
	}
	r21 = calculate.r01(n, d)[2]
	num =  num + ( d[3] - r21 * n[3] )
	den = den + n[3] * r21 * (1 - r21)
	cur[tt]  = num/den 
	}
	zscore = num/sqrt(den)
	p.value = dnorm(zscore)
	browser()
	list(zscore=zscore, p.value=p.value) 
}


interactionTest.chisq <- function(n, events, quadrant) {
# events: 0 - death; 1 - censored
n  = n+ 1
d.exp = d = rep(0,4)
sur = n  
r01.old = r01 = calculate.r01(n, d)
o1 = o3 = 0 #  d is O
e1 = e3 = 0
cur = rep(0, (length(events) -1))
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC403858/pdf/bmj32801073.pdf
for (tt in seq(length(events) -1)) {
	event = events[tt]
	quadrant.tt = quadrant[tt]
	sur[quadrant.tt] = sur[quadrant.tt] - 1
	if(event){ # censored 
		n[quadrant.tt]  = n[quadrant.tt] - 1   
	} else{ # death
		d[quadrant.tt] = d[quadrant.tt] + 1
		if(quadrant.tt %in% c(1,3)){
			r01.old = r01 
			r01 = calculate.r01(n, d)
			r01.del = r01 - r01.old
			print(r01)
			if(sum(r01.del) > 0 ){
				ratio = r01.del[2] * sur[3]/( r01.del[1] * sur[1] + r01.del[2] * sur[3] )   
				e3  = e3 + ratio * 1 # change number of deaths
				e1  = e1 + (1 - ratio) * 1 # change number of deaths

			}else{
				browser()
			}
			if(abs((d[1] + d[3]) -( e1+e3) )>.01) browser()

		}
	}

	# o1 = o1 + d[1] 
	# o3 = o3 + d[3]
	# e1 = e1 + r01[1] * n[1]
	# e3 = e3 + r01[2] * n[3]
	# cur[tt]  = num/den 
	}
	stat1 = ( d[1] - e1)^2 / e1 
	stat3 = ( d[3] - e3)^2 / e3
	stat = stat1 + stat3 
	p.value = dchisq(stat, df=1)
	stat1 = sign(d[1] - e1) * stat1
	stat3 = sign(d[3] - e3) * stat3
	browser()
	list(stat=stat, stat1=stat1, stat3=stat3, p.value=p.value) 
}
