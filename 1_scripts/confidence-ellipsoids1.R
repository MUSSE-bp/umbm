## All of this is done for the main model, X_t. 


## I took here state 4 (or 3) as initial state to get some variability, 
## taking 1 or 2 results in all observations being 1, since in the model 
## one cannot reach state 3 or 4 from state 1 or 2. 


numrepeats <-10000 ## for the simulated CLT (at the end).
## We will get a vector of the p-values of the confidence ellipse around 
## the theoretical mean (P_i) that includes the observed proportions, 
## alternatively, a vector of significance levels at which the ellipse 
## centred at the observed proportions covers the true mean. 
 


library(diversitree)

f_CLTtest<-function(pars,v_X,alpha){
    z<-sum(v_X)
    q<-length(v_X)
    v_b<-rep(1,q)

## formula for a1,a2,a3,a4
    v_a <- c(pars[["L1"]]+pars[["q12"]],
             pars[["L2"]]+pars[["q21"]],
             pars[["L3"]]+pars[["q34"]]+pars[["m32"]],
             pars[["L4"]]+pars[["q43"]])

## mean offspring matrix
    A<-rbind(c(pars[["L1"]]-pars[["q12"]],pars[["q21"]],0,0),
             c(pars[["q12"]],pars[["L2"]]-pars[["q21"]],pars[["m32"]],0),
             c(0,0,pars[["L3"]]-pars[["q34"]]-pars[["m32"]],pars[["q43"]]),
             c(0,0,pars[["q34"]],pars[["L4"]]-pars[["q43"]]))

### Right eigenvectors
    leigenRight<-eigen(A)
    vgamma<-leigenRight$values
    mEigenRightVects<-leigenRight$vectors
    
  
##  R gives eigenvectors that are already L2 normalized.
    
    mEigenRightVects<-mEigenRightVects[,order(vgamma,decreasing=TRUE)]
    vgamma<-sort(vgamma,decreasing=TRUE)
    vgammaN1<-vgamma[-1]
    
## ================================================================================

## left eigenvectors
    leigenLeft<-eigen(t(A))
    eigvals<-leigenLeft$values
    if (!isTRUE(all.equal(sort(vgamma),sort(eigvals)))){print("Eigenvalues seem different!")}
    mEigenLeftVects<-leigenLeft$vectors

    
    mEigenLeftVects<-mEigenLeftVects[,order(eigvals,decreasing=TRUE)]
    eigvals<-sort(eigvals,decreasing=TRUE)
    
    
## get dual bases of eigenvectors.. right eigenvectors need not be changed
    
    for (i in 1:4){
      mEigenLeftVects[,i]<- mEigenLeftVects[,i]/(mEigenLeftVects[,i]%*%mEigenRightVects[,i])
    }
    
    mEigenLeftVects<-matrix(c(mEigenLeftVects[,1],mEigenLeftVects[,2],
                              mEigenLeftVects[,3],mEigenLeftVects[,4] ), nrow=4,ncol=4)
    
    
    
    
## ===================================================================

## the Xi's    
    lXi<-list(list(c(1,0,0,0),c(-1,1,0,0)),
              list(c(0,1,0,0),c(1,-1,0,0)),
              list(c(0,0,1,0),c(0,0,-1,1),c(0,1,-1,0)),
              list(c(0,0,0,1),c(0,0,1,-1)))
    
    lprobXi<-list(c(pars[["L1"]],pars[["q12"]])/v_a[1],
                  c(pars[["L2"]],pars[["q21"]])/v_a[2],
                  c(pars[["L3"]],pars[["q34"]],pars[["m32"]])/v_a[3],
                  c(pars[["L4"]],pars[["q43"]])/v_a[4])
## =========================================

## matrix B
    
    lBi<-vector("list",q)
    for (i in 1:q){
	  lBi[[i]]<-matrix(0,q,q)
	  for (j in 1:length(lXi[[i]])){
	  lBi[[i]]<-lBi[[i]]+lprobXi[[i]][j]*(lXi[[i]][[j]]%*%t(lXi[[i]][[j]]))
	}
    }
    
    v1<-mEigenRightVects[,1]
    
    mB<-((v1[1]*v_a[1])*lBi[[1]])+
        ((v1[2]*v_a[2])*lBi[[2]])+
        ((v1[3]*v_a[3])*lBi[[3]])+
        ((v1[4]*v_a[4])*lBi[[4]])
    
## =========================================
    
## Matrix Sigma_I
    mSigmaI<-matrix(0,q,q)
    for (i in 2:q){
	  for (j in 2:q){
	    mSigmaI<-mSigmaI+(c((t(mEigenLeftVects[,j])%*%mB%*%mEigenLeftVects[,i]))*
	                        (mEigenRightVects[,j]%*%t(mEigenRightVects[,i])))/(vgamma[1]-vgamma[i]-vgamma[j])    
	}
    }
## =========================================    

##  Confidence ellipsoids   
## Corr 3.16 in Janson (2004)
    
  v_Corr3_16_CLT<-"Cannot calculate as regime gamma_2 < gamma_1/2 not met."
  if (vgamma[2]<vgamma[1]/2){
    
	M1<-diag(1,q,q)-(v1%*%t(v_b))/(c(v_b%*%v1))
	M2<-diag(1,q,q)-(v_b%*%t(v1))/(c(v_b%*%v1))
	Sigma_b<-M1%*%mSigmaI%*%M2/(c(v_b%*%v1))
	v_prop_mean<-v1/(c(v_b%*%v1))
	vCLT_vect<-(v_X/z-v_prop_mean)*sqrt(z)
	
	## remove the first coordinate, as all are proportions, so one dimension is defined
	type_drop<-1
	vCLT_vect<-vCLT_vect[-type_drop]
	Sigma_b<-Sigma_b[-type_drop,-type_drop,drop=FALSE]
	RSS<-vCLT_vect%*%solve(Sigma_b)%*%vCLT_vect
	p_value<-pchisq(RSS,q-1,lower.tail=FALSE)
	## Sigma_b defines the confidence ellipse for the mean(P_i) of the proportions of types 2,3,4
	## if you want another type to be dropped change the value of type_drop.
	## to get at the confidence ellipse one must look at Sigma_b's eigendecomposition
	## and use Eqs 4-7, 4-8 of Johson, Wichern, Applied Multivariate Statistical Analysis, Pearson New International Edition
	v_Corr3_16_CLT<-list(RSS=RSS,p_value=p_value,q=q,Sigma_b=Sigma_b, 
	                     eigSigma_b=eigen(Sigma_b),v_prop_mean=v_prop_mean[-type_drop],
	                     type_drop=type_drop)
    }
## =========================================
    
## Thm 3.15 in Janson (2004)
    v_Thm3_15_CLT<-NA


## =========================================

    list(Thm3_15_CLT=v_Thm3_15_CLT,Corr3_16_CLT=v_Corr3_16_CLT)
}


##############################################################################
##############################################################################
########################## astA ################################

v_astA_X<-c(108,20,7,116)
v_astA_pars<-c("L1"=26.39 ,"L2"=0.00,"L3"=229.6,"L4"=0,"q12"=176,"q21"=904.4,"q34"=407.4,"q43"=41.68,"m32"=65.53)
res<-f_CLTtest(pars=v_astA_pars,v_X=v_astA_X)
res

##############################################################################
##############################################################################
##########################   cnf1    ################################

v_cnf_X<-c(90, 38, 3, 120)
v_cnf_pars<-c("L1"= 5.840 ,"L2"=108.3,"L3"=284.8,"L4"= 8.066,"q12"=1.185,"q21"=144.6,"q34"=341.3,"q43"=14.92,"m32"=121.5)
res<-f_CLTtest(pars=v_cnf_pars,v_X=v_cnf_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q


## simulate trees
## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(m32), q34, q41(0), q42(0), q43
v_cnf_pars_diversitree<-c(5.840,108.3,284.8,8.066,
                           0,0,0,0,
                          1.185,0,0,
                          144.628,0,0,
                           0,121.5,341.3,
                           0,0,14.92)
N_cnf<-sum(v_cnf_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
## drawing root state: randomly sample 3 or 4 
## depending on proportions of strains in state 3 or 4. 
##( in principle we can also change to c(0,0,0.5,0.5) ).
## we need this due to assumption (F5).
  xpr<-c(0,0,v_cnf_X[3],v_cnf_X[4])/(sum(v_cnf_X[3:4]))
  x0<-sample(1:length(v_cnf_X),1,prob=xpr)

## generate trees  
  phy <- diversitree::tree.musse(v_cnf_pars_diversitree, N_cnf, x0=x0)
  
## we have to do this correction as sometimes one of the states is never simulated
  tip_states<-c(1,2,3,4,phy$tip.state)
  v_simX<-table(tip_states)
  v_simX<-v_simX-1
  v_RSS[i]<-NA
  
## condition on essential non-extinction:
## take only those trees in which state 3 or 4 is obtained atleast once
  
  if (sum(v_simX[3:4])>0){
    z<-sum(v_simX)
    vCLT_vect<-sqrt(z)*(v_simX[-type_drop]/z-v_prop_mean)
    v_RSS[i]<-vCLT_vect%*%solve(Sigma_b)%*%vCLT_vect
    v_p_value[i]<-pchisq(v_RSS[i],q-1,lower.tail=FALSE)
  }
}


v_p_value
v_RSS

num_removed<-length(which(is.na(v_RSS)))
num_removed
prop_num_removed<-length(which(is.na(v_RSS)))/length(v_RSS)
cnf<-v_RSS[which(!is.na(v_RSS))]
cnf
subset(cnf, cnf>10.35)
subset(cnf, cnf>11.345)
max(cnf)

## make histogram for ditribution of RSS of remaining trees
maxx<-40
maxy<-3000
h<-hist(cnf,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="cnf1",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-10.35
yy<-0
new.x<-10.35
new.y<-2000
Text<-expression(paste(D^2)==10.35)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=10.35, y=2300, label=Text, cex=0.65)



##############################################################################
##############################################################################
##########################  fimG  ################################

v_fim_X<-c(13, 115, 120, 3)
v_fim_pars<-c("L1"=4.105 ,"L2"=33.75,"L3"=96.990,"L4"=3.737,"q12"=11.743,"q21"=11.743,"q34"=3.715 ,"q43"=3.715 ,"m32"=35.184 )
res<-f_CLTtest(pars=v_fim_pars,v_X=v_fim_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(m32), q34, q41(0), q42(0), q43
v_fim_pars_diversitree<-c(4.105,33.75,96.99,3.737,
                          0,0,0,0,
                          11.743,0,0,
                          11.743,0,0,
                          0,35.184,3.715,
                          0,0,3.715)
N_fim<-sum(v_fim_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
  ## we need this due to assumption (f5)
  xpr<-c(0,0,v_fim_X[3],v_fim_X[4])/(sum(v_fim_X[3:4]))
  x0<-sample(1:length(v_fim_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_fim_pars_diversitree, N_fim, x0=x0)
  
  ## we have to do this correction as sometimes one of the states is never simulated
  tip_states<-c(1,2,3,4,phy$tip.state)
  v_simX<-table(tip_states)
  v_simX<-v_simX-1
  v_RSS[i]<-NA
  
  ## condition on essential non-extinction,
  ## take only those trees in which state 3 or 4 is obtained atleast once
  
  if (sum(v_simX[3:4])>0){
    z<-sum(v_simX)
    vCLT_vect<-sqrt(z)*(v_simX[-type_drop]/z-v_prop_mean)
    v_RSS[i]<-vCLT_vect%*%solve(Sigma_b)%*%vCLT_vect
    v_p_value[i]<-pchisq(v_RSS[i],q-1,lower.tail=FALSE)
  }
}


v_p_value
v_RSS

num_removed<-length(which(is.na(v_RSS)))
num_removed
prop_num_removed<-length(which(is.na(v_RSS)))/length(v_RSS)
fim<-v_RSS[which(!is.na(v_RSS))]
fim
subset(fim, fim>3.421)
subset(fim, fim>11.345)
max(fim)


## make histogram for ditribution of RSS of remaining trees
maxx<-40
maxy<-2000
h<-hist(fim,breaks=40, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="fimG",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-3.421
yy<-0
new.x<-3.421
new.y<-1500
Text<-expression(paste(D^2)==3.421)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=8.2, y=1700, label=Text, cex=0.65)





##############################################################################
##############################################################################
##########################  fyuA  ################################

v_fyuA_X<-c(50, 78, 74, 49)
v_fyuA_pars<-c("L1"=0.464,"L2"=51.71,"L3"=120.7,"L4"= 4.242,"q12"=0,"q21"=38.36,"q34"=47.40,"q43"=6.4445,"m32"=46.71)
res<-f_CLTtest(pars=v_fyuA_pars,v_X=v_fyuA_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(m32), q34, q41(0), q42(0), q43
v_fyuA_pars_diversitree<-c(0.464,51.71, 120.7, 4.242,
                           0,0,0,0,
                           0,0,0,
                           38.36,0,0,
                           0,46.71,47.4,
                           0,0,6.4445)
N_fyuA<-sum(v_fyuA_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
  ## we need this due to assumption F5
  xpr<-c(0,0,v_fyuA_X[3],v_fyuA_X[4])/(sum(v_fyuA_X[3:4]))
  x0<-sample(1:length(v_fyuA_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_fyuA_pars_diversitree, N_fyuA, x0=x0)
  
  ## we have to do this correction as sometimes one of the states is never simulated
  tip_states<-c(1,2,3,4,phy$tip.state)
  v_simX<-table(tip_states)
  v_simX<-v_simX-1
  v_RSS[i]<-NA
  
  ## condition on essential non-extinction,
  ## take only those trees in which state 3 or 4 is obtained atleast once
  
  if (sum(v_simX[3:4])>0){
    z<-sum(v_simX)
    vCLT_vect<-sqrt(z)*(v_simX[-type_drop]/z-v_prop_mean)
    v_RSS[i]<-vCLT_vect%*%solve(Sigma_b)%*%vCLT_vect
    v_p_value[i]<-pchisq(v_RSS[i],q-1,lower.tail=FALSE)
  }
}


v_p_value
v_RSS

num_removed<-length(which(is.na(v_RSS)))
num_removed
prop_num_removed<-length(which(is.na(v_RSS)))/length(v_RSS)
fyuA<-v_RSS[which(!is.na(v_RSS))]
fyuA
subset(fyuA, fyuA>5.655)
subset(fyuA, fyuA>11.345)
max(fyuA)


## make histogram for ditribution of RSS of remaining trees
maxx<-20
maxy<-2000
h<-hist(fyuA,breaks=28, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="fyuA",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-5.655
yy<-0
new.x<-5.655
new.y<-1300
Text<-expression(paste(D^2)==5.655)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=6, y=1500, label=Text, cex=0.65)


##############################################################################
##############################################################################
##########################  hly1  ################################

v_hly_X<-c(88, 40, 5, 118)
v_hly_pars<-c("L1"=4.630,"L2"=101.7,"L3"=265.6,"L4"=10.58,"q12"=0,"q21"=131.1,"q34"=287.9,"q43"=17.77,"m32"=110.4)
res<-f_CLTtest(pars=v_hly_pars,v_X=v_hly_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q



## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(m32), q34, q41(0), q42(0), q43
v_hly_pars_diversitree<-c(4.630,101.7,265.6,10.58,
                           0,0,0,0,
                           0,0,0,
                          131.1,0,0,
                           0,110.4,287.9,
                           0,0,17.77)
N_hly<-sum(v_hly_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
  ## we need this due to assumption (F5)
  xpr<-c(0,0,v_hly_X[3],v_hly_X[4])/(sum(v_hly_X[3:4]))
  x0<-sample(1:length(v_hly_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_hly_pars_diversitree, N_hly, x0=x0)
  
  ## we have to do this correction as sometimes one of the states is never simulated
  tip_states<-c(1,2,3,4,phy$tip.state)
  v_simX<-table(tip_states)
  v_simX<-v_simX-1
  v_RSS[i]<-NA
  
  ## condition on essential non-extinction 
  ## take only those trees in which state 3 or 4 is obtained atleast once
  
  if (sum(v_simX[3:4])>0){
    z<-sum(v_simX)
    vCLT_vect<-sqrt(z)*(v_simX[-type_drop]/z-v_prop_mean)
    v_RSS[i]<-vCLT_vect%*%solve(Sigma_b)%*%vCLT_vect
    v_p_value[i]<-pchisq(v_RSS[i],q-1,lower.tail=FALSE)
  }
}


v_p_value
v_RSS

num_removed<-length(which(is.na(v_RSS)))
num_removed
prop_num_removed<-length(which(is.na(v_RSS)))/length(v_RSS)
hly<-v_RSS[which(!is.na(v_RSS))]
hly
subset(hly, hly>8.963)
subset(hly, hly>11.345)
max(hly)

## make histogram for ditribution of RSS of remaining trees
maxx<-40
maxy<-2000
h<-hist(hly,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="hly1",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-8.963
yy<-0
new.x<-8.963
new.y<-1500
Text<-expression(paste(D^2)==8.963)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=9.6, y=1700, label=Text, cex=0.65)


##############################################################################
##############################################################################
##########################  iroN   ################################

v_iroN_X<-c(44, 84, 36, 87)
v_iroN_pars<-c("L1"=1.843,"L2"=48.50,"L3"=164.0,"L4"=0,"q12"=21.20,"q21"=60.93,"q34"=178.0,"q43"=20.55,"m32"=51.44)
res<-f_CLTtest(pars=v_iroN_pars,v_X=v_iroN_X)
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q
res


##############################################################################
##############################################################################
##########################  iutA   ################################

v_iutA_X<-c(66, 62, 77, 46)
v_iutA_pars<-c("L1"=0.000,"L2"=64.93,"L3"=126.6,"L4"=3.923,"q12"=49.07,"q21"=157.0,"q34"=43.34,"q43"=5.527,"m32"=46.21)
res<-f_CLTtest(pars=v_iutA_pars,v_X=v_iutA_X)
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q
res


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(m32), q34, q41(0), q42(0), q43
v_iutA_pars_diversitree<-c(0,64.93,126.6,3.923,
                          0,0,0,0,
                          49.07,0,0,
                          157.0,0,0,
                          0,46.21,43.34,
                          0,0,5.527)
N_iutA<-sum(v_iutA_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
  xpr<-c(0,0,v_iutA_X[3],v_iutA_X[4])/(sum(v_iutA_X[3:4]))
  x0<-sample(1:length(v_iutA_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_iutA_pars_diversitree, N_iutA, x0=x0)
  
  ## we have to do this correction as sometimes one of the states is never simulated
  tip_states<-c(1,2,3,4,phy$tip.state)
  v_simX<-table(tip_states)
  v_simX<-v_simX-1
  v_RSS[i]<-NA
  
 
  ## take only those trees in which state 3 or 4 is obtained atleast once
  
  if (sum(v_simX[3:4])>0){
    z<-sum(v_simX)
    vCLT_vect<-sqrt(z)*(v_simX[-type_drop]/z-v_prop_mean)
    v_RSS[i]<-vCLT_vect%*%solve(Sigma_b)%*%vCLT_vect
    v_p_value[i]<-pchisq(v_RSS[i],q-1,lower.tail=FALSE)
  }
}


v_p_value
v_RSS

num_removed<-length(which(is.na(v_RSS)))
num_removed
prop_num_removed<-length(which(is.na(v_RSS)))/length(v_RSS)
iutA<-v_RSS[which(!is.na(v_RSS))]
iutA
subset(iutA, iutA>5.989)
subset(iutA, iutA>11.345)
max(iutA)

## make histogram for ditribution of RSS of remaining trees
maxx<-30
maxy<-2000
h<-hist(iutA,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="iutA",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-5.989
yy<-0
new.x<-5.989
new.y<-1500
Text<-expression(paste(D^2)==5.989)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=7.3, y=1700, label=Text, cex=0.65)

##############################################################################
##############################################################################
##########################  papC  ################################

v_papC_X<-c(83, 45, 33, 90)
v_papC_pars<-c("L1"=4.314,"L2"=100.7,"L3"=164.9,"L4"=8.197,"q12"=1.915,"q21"=128.5,"q34"=152.7,"q43"=16.06,"m32"=59.41)
res<-f_CLTtest(pars=v_papC_pars,v_X=v_papC_X)
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q
res


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(m32), q34, q41(0), q42(0), q43
v_papC_pars_diversitree<-c(4.314,100.7,164.9,8.197,
                          0,0,0,0,
                          1.915,0,0,
                          128.5,0,0,
                          0,59.41,152.7,
                          0,0,16.06)

N_papC<-sum(v_papC_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
  ## we need this due to assumption (f5)
  xpr<-c(0,0,v_papC_X[3],v_papC_X[4])/(sum(v_papC_X[3:4]))
  x0<-sample(1:length(v_papC_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_papC_pars_diversitree, N_papC, x0=x0)
  
  ## we have to do this correction as sometimes one of the states is never simulated
  tip_states<-c(1,2,3,4,phy$tip.state)
  v_simX<-table(tip_states)
  v_simX<-v_simX-1
  v_RSS[i]<-NA
  
  ## 
  ## take only those trees in which state 3 or 4 is obtained atleast once
  
  if (sum(v_simX[3:4])>0){
    z<-sum(v_simX)
    vCLT_vect<-sqrt(z)*(v_simX[-type_drop]/z-v_prop_mean)
    v_RSS[i]<-vCLT_vect%*%solve(Sigma_b)%*%vCLT_vect
    v_p_value[i]<-pchisq(v_RSS[i],q-1,lower.tail=FALSE)
  }
}


v_p_value
v_RSS

num_removed<-length(which(is.na(v_RSS)))
num_removed
prop_num_removed<-length(which(is.na(v_RSS)))/length(v_RSS)
papC<-v_RSS[which(!is.na(v_RSS))]
papC
subset(papC, papC>5.989)
subset(papC, papC>11.345)
max(papC)

## make histogram for ditribution of RSS of remaining trees
maxx<-30
maxy<-2000
h<-hist(papC,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="papC",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-8.498999
yy<-0
new.x<-8.498999
new.y<-1500
Text<-expression(paste(D^2)==8.499)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=8.4, y=1700, label=Text, cex=0.65)
##############################################################################
##############################################################################
##########################  sat  ################################

v_sat_X<-c(113, 15, 46, 77)
v_sat_pars<-c("L1"=6.059,"L2"=165.6,"L3"=157.4,"L4"= 6.493,"q12"=0,"q21"=307.3,"q34"=102.4,"q43"=9.691,"m32"=65.62)
res<-f_CLTtest(pars=v_sat_pars,v_X=v_sat_X)
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q
res


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(m32), q34, q41(0), q42(0), q43
v_sat_pars_diversitree<-c(6.059,165.6,157.4,6.493,
                           0,0,0,0,
                           0,0,0,
                          307.3,0,0,
                           0,65.62,102.4,
                           0,0,9.691)
N_sat<-sum(v_sat_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).

  xpr<-c(0,0,v_sat_X[3],v_sat_X[4])/(sum(v_sat_X[3:4]))
  x0<-sample(1:length(v_sat_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_sat_pars_diversitree, N_sat, x0=x0)
  
  ## we have to do this correction as sometimes one of the states is never simulated
  tip_states<-c(1,2,3,4,phy$tip.state)
  v_simX<-table(tip_states)
  v_simX<-v_simX-1
  v_RSS[i]<-NA
  
  
  ## take only those trees in which state 3 or 4 is obtained atleast once
  
  if (sum(v_simX[3:4])>0){
    z<-sum(v_simX)
    vCLT_vect<-sqrt(z)*(v_simX[-type_drop]/z-v_prop_mean)
    v_RSS[i]<-vCLT_vect%*%solve(Sigma_b)%*%vCLT_vect
    v_p_value[i]<-pchisq(v_RSS[i],q-1,lower.tail=FALSE)
  }
}


v_p_value
v_RSS

num_removed<-length(which(is.na(v_RSS)))
num_removed
prop_num_removed<-length(which(is.na(v_RSS)))/length(v_RSS)
sat<-v_RSS[which(!is.na(v_RSS))]
sat
subset(sat, sat> 6.589)
subset(sat, sat>11.345)
max(sat)

## make histogram for ditribution of RSS of remaining trees
maxx<-40
maxy<-2000
h<-hist(sat,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="sat",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-6.588527
yy<-0
new.x<-6.588527
new.y<-1500
Text<-expression(paste(D^2)==6.589)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=9, y=1700, label=Text, cex=0.65)
