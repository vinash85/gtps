setwd('/cbcb/project2-scratch/jooslee/srescues/test')
tcga.sur = fread("Curtis_survival.txt")
tcga.sur = tcga.sur[order(V2)]
tcga.mrna = fread('Curtis_mRNA.txt')
tcga.rna = as.matrix(tcga.mrna[ ,2:1982,with=F]) 
rownames(tcga.rna) = tcga.mrna$V1
mrna.header = read.table('Curtis_mRNA.txt',nrow=1 )
colnames(tcga.rna) = as.character(unlist(mrna.header))
# tcga.rna1 = sapply(t(tcga.rna), rank)

known.sl = as.data.frame( fread('SL_list.txt', header=F))
livnat.sl = as.data.frame(fread('Livnat_SLs.txt', header=F))

# mrna.q = ncol(tcga.rna)/2 

#### divide in  quadrant
sl.pair = known.sl
for (inx in seq(nrow(sl.pairs))) {

 sl = unlist(sl.pair[inx, ])
 aa1 = tcga.rna[sl[1],] 
 aa1 = ifelse( aa1 > median(aa1) ,1,0)   
 aa2 = tcga.rna[sl[2],] 
 aa2 = ifelse( aa2 > median(aa2) ,1,0)   
 quadrant = 2*aa2 + aa1 + 1
 quadrant = quadrant[tcga.sur$V1] 
 n.temp = table(quadrant)
 n = n.temp[c("1", "2", "3", "4")]
 aa = interactionTest(n=n,events=tcga.sur$V3,quadrant=quadrant) 
 aa = interactionTest.chisq(n=n,events=tcga.sur$V3,quadrant=quadrant) 
}
 



##### run survival analysis

##solving equantion like kids
#
#
library(rSymPy)
sympy("var('ra1, ra2, rb1, rb2, r21, d12, d22, d,n11, n12, n21, n22')") # declare vars
sympy("solve([Eq(ra2+rb1-ra2*rb1, d12/n12), Eq(ra2+rb2-ra2*rb2, d22/n22), Eq(ra1+rb2-ra1*rb2, r21), Eq((n11 + n12)*rb1+ (n22 + n21) *rb2,d), Eq((n11 + n21)*ra1+ (n22 + n12) *ra2,d)],[ra1,ra2,rb1, rb2, r21])",retclass="Sym")


library(rSymPy)
sympy("var('ra1, ra2, rb1, rb2, r21, d12, d22, d,n11, n12, n21, n22')") # declare vars
sympy("solve([Eq(ra2+rb1-ra2*rb1, d12/n12), Eq(ra2+rb2-ra2*rb2, d22/n22), Eq((n11 + n12)*rb1+ (n22 + n21) *rb2,d), Eq((n11 + n21)*ra1+ (n22 + n12) *ra2,d)],[ra1,ra2,rb1, rb2, r21])",retclass="Sym")
Solve{ra2+rb1-ra2*rb1 ==  a, ra2+rb2-ra2*rb2==b, ra1+ n10*rb1+ n20 *rb2 == d, n01*ra1+ n02 * ra2 ==d, ra1,ra2,rb1, rb2}

Solve{x + y - x*y == l,
 x + z - x*z == k, 
a*y + b*z == d,
e*w + f*x == d,
w+z-w*z - u ==0
,x,y,z, u, w}

solve x+y-x y = l
x+z-x z = k
a y+b z = d
e w+f x = d
w+z-w z-u = 0  for  x, y, z, u, w

x = (a l+b k-d)/(a+b-d) and y = (-b k+b l+d (-l)+d)/(a (-l)+a-b k+b) and z = (a k-a l+d (-k)+d)/(a (-l)+a-b k+b) and 
u = (a (d (k-1)-f k l+f l-e k+e l)+(k-1) (b (d-f k)+d (-d+f+e)))/(e (a (l-1)+b (k-1))) 

and w = (a (d-f l)+b (d-f k)+d (f-d))/(e (a+b-d)) and a+b!=d and a l+b k!=a+b and b!=0


x = (a l+b k-d)/(a+b-d) and y = (b k-b l+d l-d)/(a l-a+b k-b) and z = (-a k+a l+d k-d)/(a l-a+b k-b) and
 u = (a d k-a d-a f k l+a f l-e a k+e a l+b d k-b d-b f k^2+b f k+d^2 (-k)+d^2+d f k-d f+e d k-e d)/(e (a l-a+b k-b)) 
 and w = (-a d+a f l-b d+b f k+d^2-d f)/(e (-a-b+d)) and a+b-d!=0 and a l-a+b k-b!=0 and b!=0

