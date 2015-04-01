#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>		/* constants */
#define MIN(a,b) a<b?a:b 
#define MAX(a,b) a>b?a:b 
#define ABS(a) (a>0.0) ? a:-a
// #include <iostream>
// #include <fstream>

/* Two-sample two-sided asymptotic distribution */
double
pkstwo(double x, double tol)
{
/* x[1:n] is input and output
 *
 * Compute
 *   \sum_{k=-\infty}^\infty (-1)^k e^{-2 k^2 x^2}
 *   = 1 + 2 \sum_{k=1}^\infty (-1)^k e^{-2 k^2 x^2}
 *   = \frac{\sqrt{2\pi}}{x} \sum_{k=1}^\infty \exp(-(2k-1)^2\pi^2/(8x^2))
 *
 * See e_g_ J_ Durbin (1973), Distribution Theory for Tests Based on the
 * Sample Distribution Function_  SIAM_
 *
 * The 'standard' series expansion obviously cannot be used close to 0;
 * we use the alternative series for x < 1, and a rather crude estimate
 * of the series remainder term in this case, in particular using that
 * ue^(-lu^2) \le e^(-lu^2 + u) \le e^(-(l-1)u^2 - u^2+u) \le e^(-(l-1))
 * provided that u and l are >= 1_
 *
 * (But note that for reasonable tolerances, one could simply take 0 as
 * the value for x < 0_2, and use the standard expansion otherwise_)
 *
 */
 double out, old, s, w, z;
 int k, k_max;

 k_max = (int) sqrt(2 - log(tol));

 if(x < 1) {
 	z = - (M_PI_2 * M_PI_4) / (x * x);
 	w = log(x);
 	s = 0;
 	for(k = 1; k < k_max; k += 2) {
 		s += exp(k * k * z - w);
 	}
 	out = s / M_1_SQRT_2PI;
 }
 else {
 	z = -2 * x * x;
 	s = -1;
 	k = 1;
 	old = 0;
 	out = 1;
 	while(fabs(old - out) > tol) {
 		old = out;
 		out += 2 * s * exp(z * k * k);
 		s *= -1;
 		k++;
 	}
 }
 return out;
}


using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec calculate_rx1(arma::uvec n, arma::uvec d){
	int d01 = d(0) + d(2);
	double r11, r21;
	if(d(1)==0){ //e1 = 0
		r11 = 0;
		r21 = d01 / (double) n(1); 
	}else if(d(3)==0){ // e3=0
		r11 = d01 / (double) n(0); 
		r21 = 0;
	} else{
		double temp = d01/(double) (n(3) * n(0) * d(1) + n(2) * n(1) * d(3) );
		r11 = n(3)* d(1) * temp;
		r21 = n(1)* d(3) * temp;
	}
	arma::vec r01(2);
	r01(0) = r11;
	r01(1) = r21;

	return r01;
}

using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec interactionTest_ks(arma::uvec n, arma::uvec events, arma::uvec quadrant) {
// events: 0 - death; 1 - censored
	quadrant = quadrant -1; // adjusted due to zero index 
	uvec  d(4);
	d.zeros();
	uvec sur = n;  
	vec r01(2), r01_del(2);
	r01 = calculate_rx1(n, d);
	vec r01Old = r01;
double e1 = 0, e3 = 0;
int cnt=0;
int eventLen = events.size();
mat cur(eventLen, 3);
cur.zeros();
double ratio, deno;
int quadrant_tt, event;
// # http://www_ncbi_nlm_nih_gov/pmc/articles/PMC403858/pdf/bmj32801073_pdf
for (int tt=0; tt <= eventLen-1; tt++ ) {
	event = events(tt);
	quadrant_tt = quadrant(tt);
	sur(quadrant_tt) = sur(quadrant_tt) - 1;
	if(event){ // censored 
		n(quadrant_tt)  = n(quadrant_tt) - 1;   
	} else{ // death
		d(quadrant_tt) = d(quadrant_tt) + 1;
		if((quadrant_tt ==0) | (quadrant_tt ==2)) {
			r01Old = r01; 
			r01 = calculate_rx1(n, d);
			r01_del = r01 - r01Old;
			r01_del(0) = r01_del(0)< 0? 0: r01_del(0);
			r01_del(1) = r01_del(1)< 0? 0: r01_del(1);
			// print(r01) 
			  
			if(sum(r01_del) > 0 ){
				deno = r01_del(0) * sur(0) + r01_del(1) * sur(2);
				if(deno ==0) break;
				ratio = r01_del(1) * sur(2)/(double)( deno); 
			}else{
				ratio = sur(2)/(double) ( sur(0) + sur(2) );   
			}
				e3  = e3 + ratio * 1;// # change number of deaths
				e1  = e1 + (1 - ratio) * 1; //# change number of deaths
				cur(cnt,0) = d(2); 
				cur(cnt,1) = e3; 
				cur(cnt,2) = r01(1) * n(2);
				cnt = cnt+1;


				if(abs((d(0) + d(2)) -( e1+e3) )>.01) {
					cout << "expected number of deaths not equal to observed" << endl;
					exit(-1);
				}

			}
		}

	}
	// cur = cur[1:(cnt-1), ]
	cnt++;
	cur.resize(cnt -1,3);
	// std::ofstream myfile;
  // myfile.open ("example.txt");
	// myfile << cur << endl;
	// cout << cur << endl;
  // myfile.close();
	// calculate statistics converted from ks.test.R  
	// check R-devel/src/library/stats/src/ks.test.R
	// TODO change this to account for different n_x and n_y
	int n_x = cnt-1, n_y = cnt -1;            // to avoid integer overflow
	double n_ks = n_x * n_y / (double)(n_x + n_y);
        // w = c(x, y)
	vec w(2*n_x);
	w(span(0, n_x -1)) = cur.col(0) ;
	w(span(n_x, 2*n_x-1)) = cur.col(1);
	uvec w_inx = sort_index(w);
	vec z(2*n_x); z.fill(1.0/(double)n_x);
        // uvec temp_inx = find(w_inx >= n_x);
        // z(temp_inx) = ones<vec>
	for (int ii = 0; ii < 2*n_x -1; ii++)
	{
		if(w_inx(ii) >=n_x)   z(ii) = -1.0/double(n_y);
	}
	vec z1 = cumsum(z);
        // z = cumsum(ifelse(order(w) <= n_x, 1 / n_x, - 1 / n_y))
	vec w1 = sort(w);
	int count = 0;
	for (int ii = 0; ii < 2*n_x -2; ii++)
	{
		if(w1(ii)!=w1(ii +1)){
			z(count) = z1(ii);
			count++;
		}

	}
	z(count) = z1(2*n_x -1);
	z.resize(count +1);

	// std::ofstream myfile;
	// myfile.open ("example.txt");
	// myfile << z << endl;
	// // cout << cur << endl;
	// myfile.close();
	double two_sided_stat = 0,
	greater_stat = arma::max(z),
	less_stat = - arma::min(z);
	two_sided_stat = MAX(greater_stat, less_stat);
	double tol = 1e-6;
	double two_sided_pval = 0,
	greater_pval = exp(-2 * n_ks *greater_stat*greater_stat) ,
	less_pval = exp(-2 * n_ks *less_stat*less_stat);
	two_sided_pval = 1.0 - pkstwo(sqrt(n_ks) * two_sided_stat, tol );
	vec out(6);
	out(0) = greater_stat;
	out(1) = less_stat;
	out(2) = two_sided_stat;
	out(3) = greater_pval;
	out(4) = less_pval;
	out(5) = two_sided_pval;
	return out;

}

// using namespace arma;		// shorthand
// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// arma::mat interactionTestks_gene(arma::umat gene2Mat, arma::uvec gene1) {
// 	int numgenes = gene2Mat.n_rows;
// 	mat interaction_info(numgenes, 6);
// 	uvec aa2;
// 	// interaction_info.fill(nan);
// 	for (int gene2 = 0; gene2 <=numgenes; gene2++ ) {
// 			aa2 = gene2Mat.row(gene2);
// 			cout << aa2 << endl;
// 			// inx = (!is_na(aa2)) & (!is_na(gene1))
// 			// if(sum(inx) < num_samples){
// 			// 	gene1_q = gene1[inx]
// 			// 	aa2_q = aa2[inx]
// 			// 	tcga_sur_q = tcga_sur[inx] 
// 			// 	quadrant = 2*aa2_q + gene1_q + 1
// 			// 	quadrant = quadrant[tcga_sur_q$V1] 
// 			// 	n_temp = table(quadrant)
// 			// 	n = n_temp[c("1", "2", "3", "4")]
// 			// 	interaction_info[gene2,]= unlist(interactionTest_ks(n=n,events=tcga_sur_q$V3,quadrant=quadrant)) 

// 			// }else{
// 			// 	quadrant = 2*aa2 + gene1 + 1
// 			// 	quadrant = quadrant[tcga_sur$V1] 
// 			// 	n_temp = table(quadrant)
// 			// 	n = n_temp[c("1", "2", "3", "4")]
// 			// 	interaction_info[gene2,]= unlist(interactionTest_ks(n=n,events=tcga_sur$V3,quadrant=quadrant))

// 			// }
// 	}
// 	return interaction_info;
// }

// cppFunction('(double add(NumericVector x)' {
//   int sum = x + y + z;
//   return sum;
// }')


