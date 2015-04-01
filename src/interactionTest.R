calculate.rx1 <- function(n, d){
	d01 = d[1] + d[3]

	if(d[2]==0){ #e1 = 0
		r11 = 0
		r21 = d01 / n[3] 
	}else if(d[4]==0){ # e3=0
		r11 = d01 / n[1] 
		r21 = 0
	} else{
		temp = d01/(n[4] * n[1] * d[2] + n[3] * n[2] * d[4] )
		r11 = n[4]* d[2] * temp
		r21 = n[2]* d[4] * temp
	}
	# print(c(r11, r21))
  return(c(r11, r21))
  # return(c(r21, r11))
}
calculate.r01 <- function(n, d){
	N = sum(n) 
	D = sum(d)
	  a = n[1] + n[2]
	  b = n[3] + n[4]
	  g = n[1] + n[3]
	  f = n[2] + n[4]
	  l = d[2]/n[2]
	  k = d[4]/n[4]
	x=(a*l+b*k-D)/(a+b-D)
	y=(-b*k+b*l+D*(-l)+D)/(a*(-l)+a-b*k+b)
	z=(-a*k+a*l+D*k-D)/(a*l-a+b*k-b)
	w=(a*D-a*f*l+b*D-b*f*k-D^2+D*f)/(g*(a+b-D))
	ra2 = x
	rb1 = y
	rb2 = z
	ra1 = w
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
	p.value = 1- pnorm(zscore)
	browser()
	list(zscore=zscore, p.value=p.value) 
}


interactionTest.ks <- function(n, events, quadrant, plot=F) {
# events: 0 - death; 1 - censored
# n  = n+ 1
d.exp = d = rep(0,4)
sur = n  
r01.old = r01 = calculate.rx1(n, d)
o1 = o3 = 0 #  d is O
e1 = e3 = 0
cnt=1
cur = matrix(0, nrow=(length(events) -1), ncol=3)
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
			r01 = calculate.rx1(n, d)
			r01.del = r01 - r01.old
			r01.del[r01.del < 0] =0
			# print(r01) 
			if(sum(r01.del) > 0 ){
				ratio = r01.del[2] * sur[3]/( r01.del[1] * sur[1] + r01.del[2] * sur[3] )   
			}else{
				ratio = sur[3]/( sur[1] + sur[3] )   
				# print("r01 is equal to zero")
				# browser()
			}
				e3  = e3 + ratio * 1 # change number of deaths
				e1  = e1 + (1 - ratio) * 1 # change number of deaths
				cur[cnt,] = c(d[3], e3, r01[2] * n[3])
				cnt = cnt+1


			if(abs((d[1] + d[3]) -( e1+e3) )>.01) {
				print("expected number of deaths not equal to observed")
			
				exit()
			}

		}
	}

	}
	cur = cur [1:(cnt-1), ]
	stat1 = ( d[1] - e1)*(d[1] - e1) / e1 
	stat3 = ( d[3] - e3) *(d[3] - e3)/ e3
	stat = stat1 + stat3 
	p.value = 1- pchisq(stat, df=1)
	stat1 = sign(d[1] - e1) * stat1
	stat3 = sign(d[3] - e3) * stat3
	two.sided.p = (ks.test(cur[,1], cur[,2], alternative="two.sided"))

	greater.p = (ks.test(cur[,1], cur[,2], alternative="greater"))
	less.p = (ks.test(cur[,1], cur[,2], alternative="less"))	# browser()
	# list(stat=stat, stat1=stat1, stat3=stat3, p.value=p.value, cur=cur, greater.p=greater.p, less.p=less.p, two.sided.p=two.sided.p)
	if(plot){
		out = list(
			less.p=less.p$p.value, 
			greater.p=greater.p$p.value, 
			two.sided.p=two.sided.p$p.value, 
			less.stat=less.p$statistic, 
			greater.stat=greater.p$statistic, 
			two.sided.stat=two.sided.p$statistic,
			cur = cur,
			d = d

			)
	} else {
		out = list(
			less.p=less.p$p.value, 
			greater.p=greater.p$p.value, 
			two.sided.p=two.sided.p$p.value, 
			less.stat=less.p$statistic, 
			greater.stat=greater.p$statistic, 
			two.sided.stat=two.sided.p$statistic 
			)
	}
	return(out)
}


plotInteractionTest.ks = function(cur) {
	require(ggplot2)
	n = nrow(cur)

	dt = rbind(
		data.table(fraction=cur[,1], inx = seq(n), exp = "Observed"),
		data.table(fraction=cur[,2],inx = seq(n), exp = "expected"),
		data.table(fraction=cur[,3],inx = seq(n), exp = "running-mean")
		)
	p  = ggplot(data=dt, aes(y=fraction, x = inx, color=exp)) + geom_line()
	p
}









