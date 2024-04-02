#' Fine-Mapping and Estimation of Causal Variants and One-Step Infinitesimal Effect Estimation Using SuSiE
#'
#' This function extends the SuSiE (Sum of Single Effects) method for fine-mapping to not only identify causal variants with significant effect sizes but also to perform a one-step estimation of the infinitesimal effect based on the residuals of the fine-mapping results. It employs REML (Restricted Maximum Likelihood) for the estimation of infinitesimal effects and for determining the optimal variance of these effects, following the identification of causal variants. This one-step approach aims to streamline the process by directly estimating the infinitesimal effect after the causal variants have been identified, enhancing the efficiency and accuracy of genetic effect estimation.
#'
#' @param z Z scores of GWAS effect sizes.
#' @param R LD (Linkage Disequilibrium) matrix of variants.
#' @param n Sample size of GWAS data.
#' @param L Number of single effects to consider, default is 10.
#' @param pip.threshold A threshold to determine which variants are considered causal, default is 0.25.
#' @param iter Maximum number of iterations for the one-step estimation process, default is 15.
#' @param inner.iter Number of iterations for the inner loop of the one-step estimation process, default is 5.
#' @param score.test Perform score test of variance component or not, default to F.
#' @param pleiotropy.keep The indices of variants to input for fine-mapping, default to "ALL".
#' @param estimate_residual_variance Boolean flag to indicate whether to estimate the residual variance in the one-step process, default is FALSE.
#' @return A list containing the following elements:
#'   - \code{eta}: Linear predictor, which is the sum of the causal effect (\code{beta}) and the infinitesimal effect (\code{alpha}). It represents the total genetic effect of variants on the trait, combining the direct effects of causal variants and the background genetic effects.
#'   - \code{beta}: Causal effect. These are the effect estimates for variants identified with significant effect sizes through the SuSiE method, representing the direct impact of these variants on the trait.
#'   - \code{alpha}: Infinitesimal effect. This represents the aggregate effect of all other genetic variants in the background, excluding the identified causal variants. It can be considered as the unexplained genetic variability in the model.
#'   - \code{var.inf}: Variance of the infinitesimal effect. Estimated via REML in a one-step process, this reflects the degree of variation in the infinitesimal effect across the genetic variation.
#'   - \code{pip}: Posterior Inclusion Probabilities for the estimates of causal effects. These are the probabilities calculated by the SuSiE method for each variant, assessing the likelihood that a variant is causal.
#'   - \code{fit.susie}: Output from the SuSiE method. This is the complete result object from the SuSiE method after performing fine-mapping and the one-step infinitesimal effect estimation, containing detailed information about the selected model and estimated parameters, allowing for further analysis and interpretation of the results.
#' @import CppMatrix
#' @importFrom susieR.nowarning susie_rss
#' @importFrom Matrix bdiag
#' @examples
#' # Example usage:
#' # Assuming z, R, and n are defined:
#' result <- susie_ios_onestep(z, R, n, L, pip.threshold, iter, inner.iter, estimate_residual_variance)
#' @export
#'
#' @details
#' The function begins with the `susie_rss` method for fine-mapping based on the input Z scores, LD matrix, and sample size to identify causal variants. After identifying these variants, it utilizes a one-step process to estimate the infinitesimal effect based on the residuals. This process involves adjusting the estimation of the infinitesimal effect to account for the identified causal variants, thereby refining the overall genetic effect estimation. The inclusion of a one-step estimation process for the infinitesimal effect aims to enhance the precision and efficiency of genetic architecture estimation by directly incorporating the fine-mapping results into the infinitesimal effect estimation.
susie_ios_onestep <- function(z, R, n, L = 10, pip.threshold = 0.25, iter = 5, inner.iter = 3, score.test = F, pleiotropy.keep = "ALL", estimate_residual_variance = F) {

if(pleiotropy.keep[1]=="ALL"){
pleiotropy.keep=c(1:length(z))
}

var.inf=0.5
alpha=beta=z*0
D=diag(R)
fiteigen=matrixEigen(R)
U=fiteigen$vector
Gamma=fiteigen$values
m=length(z)
pip=beta*0
fit=susie_rss(z=z[pleiotropy.keep],R=R[pleiotropy.keep,pleiotropy.keep],n=n,L=max(1,L),estimate_residual_variance=estimate_residual_variance)
beta[pleiotropy.keep]=coef(fit)[-1]*(fit$pip>=pip.threshold)*sqrt(n)*(L!=0)
indvalid=which(beta==0)
indpleiotropy=which(beta!=0)
pip[pleiotropy.keep]=fit$pip

if(length(indpleiotropy)>0){
G=bdiag(R,R[indpleiotropy,indpleiotropy])
G=as.matrix(G)
G[1:m,-c(1:m)]=R[,indpleiotropy]
G[-c(1:m),1:m]=R[indpleiotropy,]
g=G[1,]*0
g[1:m]=1
for(i in 1:iter){
Hinv=matrixGeneralizedInverse(G+1/var.inf*diag(g))
df=sum(diag(Hinv)[1:m])
zeta=matrixVectorMultiply(Hinv,c(z,z[indpleiotropy]))
alpha=zeta[1:m]
for(j in 1:inner.iter){
Hinv=matrixGeneralizedInverse(G+1/var.inf*diag(g))
df=sum(diag(Hinv)[1:m])
var.inf=(sum(alpha^2)+df)/m
}
}
}

if(length(indpleiotropy)==0){
G=LD
Hinv=matrixMultiply(U,t(U)*(1/(Gamma+1/var.inf)))
for(i in 1:iter){
alpha=c(matrixVectorMultiply(Hinv,z))
for(j in 1:inner.iter){
Hinv=matrixMultiply(U,t(U)*(1/(Gamma+1/var.inf)))
df=sum(diag(Hinv))
var.inf=(sum(alpha^2)+df)/m
}
}
}

res=c(z-matrixVectorMultiply(R,beta))
if(score.test==T){
pv=inf_test(res.inf=res,LD=R,Theta=matrixInverse(R),A=R[,which(beta!=0)])
}else{
pv=1
}
return(list(eta=alpha+beta,beta=beta,alpha=alpha,df.eta=sum(diag(matrixMultiply(Hinv,G))),var.inf=var.inf,pv=pv,pip=pip,fit.susie=fit))
}
