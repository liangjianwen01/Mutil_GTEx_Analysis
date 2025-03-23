confint<-function(x,sigma=-1,alpha=0.05)
{
  n<-length(x)
  xb<-mean(x)
  if(sigma>=0)
  {
    tmp<-sigma/sqrt(n)*qnorm(1-alpha/2);df<-n
  }
  else{
    tmp<-sd(x)/sqrt(n)*qt(1-alpha/2,n-1);df<- n-1
  }
  data.frame(mean=xb,df=df,a=xb-tmp,b=xb+tmp)
}

Adrenal_tissues <- read.csv('~/Rpro/Mutil-GTEx/Adrenal_tissues.csv')[,-1]
qnorm(0.025,mean = mean(Adrenal_tissues$mR), sd(Adrenal_tissues$mR))#177
qnorm(0.975,mean = mean(Adrenal_tissues$mR), sd(Adrenal_tissues$mR))#198
qnorm(0.025,mean = mean(Adrenal_tissues$mG), sd(Adrenal_tissues$mG))#109
qnorm(0.975,mean = mean(Adrenal_tissues$mG), sd(Adrenal_tissues$mG))#147
qnorm(0.025,mean = mean(Adrenal_tissues$mB), sd(Adrenal_tissues$mB))#148
qnorm(0.975,mean = mean(Adrenal_tissues$mB), sd(Adrenal_tissues$mB))#179

Pituitary_tissues <- read.csv('~/Rpro/Mutil-GTEx/Pituitary_tissues.csv')
qnorm(0.025,mean = mean(Pituitary_tissues$mR), sd(Pituitary_tissues$mR))#156
qnorm(0.975,mean = mean(Pituitary_tissues$mR), sd(Pituitary_tissues$mR))#182
qnorm(0.025,mean = mean(Pituitary_tissues$mG), sd(Pituitary_tissues$mG))#66
qnorm(0.975,mean = mean(Pituitary_tissues$mG), sd(Pituitary_tissues$mG))#124
qnorm(0.025,mean = mean(Pituitary_tissues$mB), sd(Pituitary_tissues$mB))#106
qnorm(0.975,mean = mean(Pituitary_tissues$mB), sd(Pituitary_tissues$mB))#158

appendage_tissues <- read.csv('~/Rpro/Mutil-GTEx/appendage_tissues.csv')
qnorm(0.025,mean = mean(appendage_tissues$mR), sd(appendage_tissues$mR))#187
qnorm(0.975,mean = mean(appendage_tissues$mR), sd(appendage_tissues$mR))#205
qnorm(0.025,mean = mean(appendage_tissues$mG), sd(appendage_tissues$mG))#93
qnorm(0.975,mean = mean(appendage_tissues$mG), sd(appendage_tissues$mG))#135
qnorm(0.025,mean = mean(appendage_tissues$mB), sd(appendage_tissues$mB))#133
qnorm(0.975,mean = mean(appendage_tissues$mB), sd(appendage_tissues$mB))#166

gj_tissues <- read.csv('~/Rpro/Mutil-GTEx/gj_tissues.csv')
qnorm(0.025,mean = mean(gj_tissues$mR), sd(gj_tissues$mR))#203
qnorm(0.975,mean = mean(gj_tissues$mR), sd(gj_tissues$mR))#224
qnorm(0.025,mean = mean(gj_tissues$mG), sd(gj_tissues$mG))#104
qnorm(0.975,mean = mean(gj_tissues$mG), sd(gj_tissues$mG))#170
qnorm(0.025,mean = mean(gj_tissues$mB), sd(gj_tissues$mB))#140
qnorm(0.975,mean = mean(gj_tissues$mB), sd(gj_tissues$mB))#191

sc_tissues <- read.csv('~/Rpro/Mutil-GTEx/sigmoidcolon_tissues.csv')
qnorm(0.025,mean = mean(sc_tissues$mR), sd(sc_tissues$mR))#182
qnorm(0.975,mean = mean(sc_tissues$mR), sd(sc_tissues$mR))#214
qnorm(0.025,mean = mean(sc_tissues$mG), sd(sc_tissues$mG))#51
qnorm(0.975,mean = mean(sc_tissues$mG), sd(sc_tissues$mG))#145
qnorm(0.025,mean = mean(sc_tissues$mB), sd(sc_tissues$mB))#90
qnorm(0.975,mean = mean(sc_tissues$mB), sd(sc_tissues$mB))#169
