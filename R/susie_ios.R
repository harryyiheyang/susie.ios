#' Fine-Mapping and Estimation of Causal Variants Using SuSiE
#'
#' This function utilizes the SuSiE (Sum of Single Effects) method for fine-mapping to identify causal variants with significant effect sizes. It employs REML (Restricted Maximum Likelihood) for the estimation of infinitesimal effects and for determining the optimal variance of these effects. Additionally, the function utilizes a score test to assess whether the variance of the infinitesimal effect is zero. The function has been updated to include parameters for controlling the precision of the estimation process and the option to estimate residual variance.
#'
#' @param z Z scores of GWAS effect sizes.
#' @param R LD (Linkage Disequilibrium) matrix of variants.
#' @param n Sample size of GWAS data.
#' @param L Number of single effects to consider, default is 10.
#' @param pip.threshold A threshold to determine which variants are considered causal, default is 0.5.
#' @param max.iter Maximum number of iterations for estimating the infinitesimal effect, default is 15.
#' @param max.eps The maximum epsilon for convergence in the iterative process, default is 0.001.
#' @param inner.iter Number of iterations for the inner loop of the estimation process, default is 5.
#' @param score.test Perform score test of variance component or not, default to F.
#' @param pleiotropy.keep The indices of variants to input for fine-mapping, default to "ALL".
#' @param estimate_residual_variance Boolean flag to indicate whether to estimate the residual variance, default is FALSE.
#' @return A list containing the following elements:
#'   - \code{eta}: Linear predictor, which is the sum of the causal effect (\code{beta}) and the infinitesimal effect (\code{alpha}). It represents the total genetic effect of variants on the trait, combining the direct effects of causal variants and the background genetic effects.
#'   - \code{beta}: Causal effect. These are the effect estimates for variants identified with significant effect sizes through the SuSiE method, representing the direct impact of these variants on the trait.
#'   - \code{alpha}: Infinitesimal effect. This represents the aggregate effect of all other genetic variants in the background, excluding the identified causal variants. It can be considered as the unexplained genetic variability in the model.
#'   - \code{var.inf}: Variance of the infinitesimal effect. Estimated via REML, this reflects the degree of variation in the infinitesimal effect across the genetic variation.
#'   - \code{pv}: P-value from the score test. This value is used to assess whether the variance of the infinitesimal effect is significantly non-zero, serving as an important metric to test if the model adequately explains the genetic variability.
#'   - \code{pip}: Posterior Inclusion Probabilities for the estimates of causal effects. These are the probabilities calculated by the SuSiE method for each variant, assessing the likelihood that a variant is causal.
#'   - \code{fit.susie}: Output from the SuSiE method. This is the complete result object from the SuSiE method after performing fine-mapping, containing detailed information about the selected model and estimated parameters, allowing for further analysis and interpretation of the results.
#' @import CppMatrix
#' @importFrom susieR.nowarning susie_rss
#' @importFrom Matrix bdiag
#' @examples
#' # Example usage:
#' # Assuming z, R, and n are defined:
#' result <- susie_ios(z, R, n, L, pip.threshold, max.iter, max.eps, inner.iter, estimate_residual_variance)
#' @export
#'
#' @details
#' The function starts by using the `susie_rss` function for fine-mapping based on the input Z scores, LD matrix, and sample size to identify causal variants. It then filters these variants based on the PIP threshold to calculate their effect sizes. For the estimation of infinitesimal effects, if the score test's p-value is below the specified threshold, an iterative process is used to estimate the optimal variance of these effects. If the p-value does not meet the threshold, the variance of the infinitesimal effect is set to zero. This approach aims to enhance the precision of fine-mapping by selecting causal variants based on their effect sizes and accurately estimating the underlying genetic architecture through REML. The inclusion of parameters for controlling the precision of the estimation process and the option to estimate residual variance allows for more flexibility and accuracy in the analysis.
susie_ios <- function(z, R, n, L = 10, pip.threshold = 0.5, max.iter = 10, max.eps = 0.001, inner.iter = 5, score.test = F, pleiotropy.keep = "ALL" , estimate_residual_variance = F) {

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

fit=susie_rss(z=z[pleiotropy.keep],R=R[pleiotropy.keep,pleiotropy.keep],n=n,L=L,estimate_residual_variance=estimate_residual_variance)
beta[pleiotropy.keep]=coef(fit)[-1]*(fit$pip>=pip.threshold)*sqrt(n)
beta1=beta*0
error=1
iter=0

while(error>max.eps&iter<max.iter){
beta1=beta
beta=beta*0
s=c(z-matrixVectorMultiply(R,alpha))
fit=susie_rss(z=s[pleiotropy.keep],R=R[pleiotropy.keep,pleiotropy.keep],n=n,L=L,estimate_residual_variance=estimate_residual_variance)
beta[pleiotropy.keep]=coef(fit)[-1]*(fit$pip>=pip.threshold)*sqrt(n)
indvalid=which(beta==0)
indpleiotropy=which(beta!=0)

for(j in 1:inner.iter){
if(length(indpleiotropy)>0){
G=bdiag(R,R[indpleiotropy,indpleiotropy])
G=as.matrix(G)
G[1:m,-c(1:m)]=R[,indpleiotropy]
G[-c(1:m),1:m]=R[indpleiotropy,]
g=G[1,]*0
g[1:m]=1
Hinv=matrixGeneralizedInverse(G+1/var.inf*diag(g))
df=sum(diag(Hinv)[1:m])
zeta=matrixVectorMultiply(Hinv,c(z,z[indpleiotropy]))
alpha=zeta[1:m]
var.inf=(sum(alpha^2)+df)/m
}else{
Hinv=matrixMultiply(U,t(U)*(1/(Gamma+1/var.inf)))
alpha=c(matrixVectorMultiply(Hinv,z))
df=sum(diag(Hinv))
var.inf=(sum(alpha^2)+df)/m
}
}
if(iter>3){error=max(abs(beta-beta1))}
iter=iter+1
}

res=c(z-matrixVectorMultiply(R,beta))
if(score.test==T){
pv=inf_test(res.inf=res,LD=R,Theta=matrixInverse(R),A=R[,which(beta!=0)])
}else{
pv=1
}
return(list(eta=alpha+beta,beta=beta,alpha=alpha,var.inf=var.inf,pv=pv,pip=fit$pip,fit.susie=fit))
}
