library(data.table)
library(glue)
library(CppMatrix)
library(susieR)
library(Matrix)
n=1e6
m=200
R=kronecker(matrix(0.1,m/4,m/4)+diag(m/4)*0.9,toeplitz(c(1,0.5,0.25,0.125)))
RC=chol(R)
A=matrix(0,100,4)
B=matrix(0,100,3)
C=D=matrix(0,100,m)
temp_dir="~/susie_inf"
write.table(R,glue("{temp_dir}/LD.txt"),row.names=F,col.names=F,quote=F)
system(glue("gzip {temp_dir}/LD.txt"))

for(i in 1:100){
  beta=rbinom(m,1,0.025)
  beta=beta*sqrt(0.00025/sum(beta^2))
  alpha=rnorm(m,0,1)
  alpha=alpha*sqrt(0.00025/sum(alpha^2))
  C[i,]=beta
  D[i,]=alpha
  z=c(R%*%(beta+alpha)*sqrt(n)+RC%*%rnorm(m,0,1))
  write.table(data.frame(z=z),glue("{temp_dir}/zscore_{i}.txt"),row.names=F,quote=F)
  fit1=susieR::susie_rss(z=z,R=R,n=n,L=5)
  fit2=susie_ios(z=z,R=R,n=n,L=5)
  fit3=susie_ios_onestep(z=z,R=R,n=n,L=5)
  A[i,1:3]=c(cor((fit1$pip>0.5)*coef(fit1)[-1]*sqrt(n),beta),cor(fit2$beta,beta),cor(fit3$beta,beta))
  B[i,1:2]=c(cor(fit2$alpha,alpha),cor(fit3$alpha,alpha))
}

for(i in 1:100){
  w=readr::read_tsv(glue("~/susie_inf/result_{i}.susieinf.bgz"))
  A[i,4]=cor(w$post_mean,C[i,])
  B[i,3]=cor(w$alpha,D[i,])
}
