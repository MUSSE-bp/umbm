##third model, Y_t. 

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
           pars[["L4"]]+pars[["q43"]]+pars[["m42"]])
  
  ## mean offspring matrix
  A<-rbind(c(pars[["L1"]]-pars[["q12"]],pars[["q21"]],0,0),
           c(pars[["q12"]],pars[["L2"]]-pars[["q21"]],pars[["m32"]],pars[["m42"]]),
           c(0,0,pars[["L3"]]-pars[["q34"]]-pars[["m32"]],pars[["q43"]]),
           c(0,0,pars[["q34"]],pars[["L4"]]-pars[["q43"]]-pars[["m42"]]))
  
  ### Right eigenvectors
  leigenRight<-eigen(A)
  vgamma<-leigenRight$values
  mEigenRightVects<-leigenRight$vectors
  
  
  
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
            list(c(0,0,0,1),c(0,0,1,-1),c(0,1,0,-1)))
  
  lprobXi<-list(c(pars[["L1"]],pars[["q12"]])/v_a[1],
                c(pars[["L2"]],pars[["q21"]])/v_a[2],
                c(pars[["L3"]],pars[["q34"]],pars[["m32"]])/v_a[3],
                c(pars[["L4"]],pars[["q43"]],pars[["m42"]])/v_a[4])
  ## =========================================
  
  ## ## matrix B
  
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
    
    ## remove the first coordinate, as all are proportions, so one dimension is defined,
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
v_astA_pars<-c("L1"=3.031,"L2"=138.8,"L3"=0,"L4"=86.05,"q12"=27.86,"q21"=448.9,"q34"=8.210,"q43"=18.47,"m32"=0,"m42"=27.88)
res<-f_CLTtest(pars=v_astA_pars,v_X=v_astA_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q


## simulate trees
## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32, q34, q41(0), m42, q43
v_astA_pars_diversitree<-c(3.031,138.8,0,86.05,
                           0,0,0,0,
                           27.86,0,0,
                           448.9,0,0,
                           0,0,8.210,
                           0,27.88,18.47)
N_astA<-sum(v_astA_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
  ## we need this due to assumption (F5) 
  xpr<-c(0,0,v_astA_X[3],v_astA_X[4])/(sum(v_astA_X[3:4]))
  x0<-sample(1:length(v_astA_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_astA_pars_diversitree, N_astA, x0=x0)
  
  ## we have to do this correction as sometimes one of the states is never simulated
  tip_states<-c(1,2,3,4,phy$tip.state)
  v_simX<-table(tip_states)
  v_simX<-v_simX-1
  v_RSS[i]<-NA
  
  ## conditione on essential non-extinction,
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
astA<-v_RSS[which(!is.na(v_RSS))]
astA #6757
subset(astA, astA>19.35) #283
subset(astA, astA>11.345) #409
max(astA) #44.14641

## make histogram for ditribution of RSS of remaining trees
maxx<-50
maxy<-3600
h<-hist(astA,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, 
     col=c("white","black")[ccat], xlim = c(0,maxx), 
     ylim = c(0,maxy), main="astA",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), 
     labels=as.character(seq(from=0,to=maxy,by=1000)), 
     cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), 
     labels=as.character(seq(from=0,to=maxx,by=10)), 
     cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-19.35
yy<-0
new.x<-19.35
new.y<-2500
Text<-expression(paste(D^2)==19.35)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=19.35, y=2800, label=Text, cex=0.65)


##############################################################################
##############################################################################
##########################   cnf1    ################################

v_cnf_X<-c(90, 38, 3, 120)
v_cnf_pars<-c("L1"= 0 ,"L2"=24.84,"L3"=3.768,"L4"= 107.5,"q12"=2.343,"q21"=136.6,"q34"=2.027,"q43"=2.706,"m32"=2.053,"m42"=33.95)
res<-f_CLTtest(pars=v_cnf_pars,v_X=v_cnf_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q


## simulate trees
## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(m32), q34, q41(0), m42, q43
v_cnf_pars_diversitree<-c(0,24.84,3.768,107.5,
                          0,0,0,0,
                          2.343,0,0,
                          136.6,0,0,
                          0,2.053,2.027,
                          0,33.95,2.706)
N_cnf<-sum(v_cnf_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
 
  xpr<-c(0,0,v_cnf_X[3],v_cnf_X[4])/(sum(v_cnf_X[3:4]))
  x0<-sample(1:length(v_cnf_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_cnf_pars_diversitree, N_cnf, x0=x0)
  
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
cnf<-v_RSS[which(!is.na(v_RSS))]
cnf #6809
subset(cnf, cnf>14.722) #31
subset(cnf, cnf>11.345) #79
max(cnf) #29.55134

## make histogram for ditribution of RSS of remaining trees
maxx<-30
maxy<-2000
h<-hist(cnf,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="cnf1",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-14.722
yy<-0
new.x<-14.722
new.y<-1500
Text<-expression(paste(D^2)==14.72)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=14.722, y=1700, label=Text, cex=0.65)



##############################################################################
##############################################################################
##########################  fimG  ################################

v_fim_X<-c(13, 115, 120, 3)
v_fim_pars<-c("L1"=4.487 ,"L2"=35.10,"L3"=103.2,"L4"=4.034,"q12"=0,"q21"=11.17 ,"q34"=1.832 ,"q43"=2.014,"m32"=30.41,"m42"=2.125 )
res<-f_CLTtest(pars=v_fim_pars,v_X=v_fim_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32, q34, q41(0), q42, q43
v_fim_pars_diversitree<-c(4.487,35.10,103.2,4.034,
                          0,0,0,0,
                          0,0,0,
                          11.17,0,0,
                          0,30.41,1.832,
                          0,2.125,2.014)
N_fim<-sum(v_fim_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
  
  xpr<-c(0,0,v_fim_X[3],v_fim_X[4])/(sum(v_fim_X[3:4]))
  x0<-sample(1:length(v_fim_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_fim_pars_diversitree, N_fim, x0=x0)
  
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
fim<-v_RSS[which(!is.na(v_RSS))]
fim #6993
subset(fim, fim>4.203041) #1276
subset(fim, fim>11.345) #221
max(fim) # 146.7622


## make histogram for ditribution of RSS of remaining trees
maxx<-150
maxy<-6800
h<-hist(fim,breaks=40, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="fimG",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=30), labels=as.character(seq(from=0,to=maxx,by=30)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-4.203041
yy<-0
new.x<-4.203041
new.y<-6400
Text<-expression(paste(D^2)==4.203)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=20, y=6700, label=Text, cex=0.65)





##############################################################################
##############################################################################
##########################  fyuA  ################################

v_fyuA_X<-c(50, 78, 74, 49)
v_fyuA_pars<-c("L1"=3.705,"L2"=57.57,"L3"=116.2,"L4"= 0,"q12"=8.254,"q21"=56.65 ,"q34"=69.75,"q43"=12.72,"m32"=36.23,"m42"=0)
res<-f_CLTtest(pars=v_fyuA_pars,v_X=v_fyuA_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q



##############################################################################
##############################################################################
##########################  hly1  ################################

v_hly_X<-c(88, 40, 5, 118)
v_hly_pars<-c("L1"=4.768,"L2"=109.4,"L3"=0,"L4"=86.64,"q12"=0.823,"q21"=145.7,"q34"=8.046,"q43"=17.18,"m32"=0,"m42"=27.30)
res<-f_CLTtest(pars=v_hly_pars,v_X=v_hly_X)
res
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q



## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32, q34, q41(0), q42, q43
v_hly_pars_diversitree<-c(4.768,109.4,0,86.64,
                          0,0,0,0,
                          0.823,0,0,
                          145.7,0,0,
                          0,0,8.046,
                          0,27.30,17.18)
N_hly<-sum(v_hly_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
 
  xpr<-c(0,0,v_hly_X[3],v_hly_X[4])/(sum(v_hly_X[3:4]))
  x0<-sample(1:length(v_hly_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_hly_pars_diversitree, N_hly, x0=x0)
  
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
hly<-v_RSS[which(!is.na(v_RSS))]
hly #6845
subset(hly, hly>20.81) #85
subset(hly, hly>11.345) #175
max(hly) #36.58013

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


xx<-20.81
yy<-0
new.x<-20.81
new.y<-1500
Text<-expression(paste(D^2)==20.81)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=20.81, y=1700, label=Text, cex=0.65)


##############################################################################
##############################################################################
##########################  iroN   ################################

v_iroN_X<-c(44, 84, 36, 87)
v_iroN_pars<-c("L1"=1.902,"L2"=48.60 ,"L3"= 164.7,"L4"=0,"q12"=20.77,"q21"=60.34,"q34"=181.2,"q43"=20.79,"m32"=51.5,"m42"=0)
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
v_iutA_pars<-c("L1"=0.000,"L2"=71.13,"L3"=155.1,"L4"=7.056,"q12"=38.11,"q21"=148.3,"q34"=47.56,"q43"=1.058,"m32"=9.561,"m42"=9.561)
res<-f_CLTtest(pars=v_iutA_pars,v_X=v_iutA_X)
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q
res


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32, q34, q41(0), q42, q43
v_iutA_pars_diversitree<-c(0, 71.13, 155.1, 7.056,
                           0,0,0,0,
                           38.11,0,0,
                           148.3,0,0,
                           0,9.561,47.56,
                           0,9.561,1.058)
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
iutA #5421
subset(iutA, iutA>93.70) #701
subset(iutA, iutA>11.345) #866
max(iutA) #791.6552

## make histogram for ditribution of RSS of remaining trees
maxx<-800
maxy<-6000
h<-hist(iutA,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="iutA",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=2000), labels=as.character(seq(from=0,to=maxy,by=2000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=200), labels=as.character(seq(from=0,to=maxx,by=200)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-93.70
yy<-0
new.x<-93.70
new.y<-3500
Text<-expression(paste(D^2)==93.70)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=135, y=4000, label=Text, cex=0.65)

##############################################################################
##############################################################################
##########################  papC  ################################

v_papC_X<-c(83, 45, 33, 90)
v_papC_pars<-c("L1"=4.493,"L2"=99.45,"L3"=166.2,"L4"=0,"q12"=2.081,"q21"=126.5,"q34"=203.0 ,"q43"=25.97,"m32"=49.39,"m42"=0.027)
res<-f_CLTtest(pars=v_papC_pars,v_X=v_papC_X)
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q
res


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(0), q34, q41(0), q42, q43
v_papC_pars_diversitree<-c(4.493,99.45,166.2,0,
                           0,0,0,0,
                           2.081,0,0,
                           126.5,0,0,
                           0,49.39,203,
                           0,0.027,25.97)

N_papC<-sum(v_papC_X)
v_p_value<-rep(0,numrepeats)
v_RSS<-rep(0,numrepeats)
l_simX<-vector("list",numrepeats)

for (i in 1:numrepeats){
  ## drawing root state: randomly sample 3 or 4 
  ## depending on proportions of strains in state 3 or 4. 
  ##( in principle we can also change to c(0,0,0.5,0.5) ).
 
  xpr<-c(0,0,v_papC_X[3],v_papC_X[4])/(sum(v_papC_X[3:4]))
  x0<-sample(1:length(v_papC_X),1,prob=xpr)
  
  ## generate trees  
  phy <- diversitree::tree.musse(v_papC_pars_diversitree, N_papC, x0=x0)
  
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
papC<-v_RSS[which(!is.na(v_RSS))]
papC # 7033
subset(papC, papC>9.619) #161
subset(papC, papC>11.345)#103
max(papC)#38.9086

## make histogram for ditribution of RSS of remaining trees
maxx<-40
maxy<-2000
h<-hist(papC,breaks=30, xlim = c(0,maxx),ylim = c(0,maxy))
ccat = cut(h$breaks, c(-Inf, 11.345, maxx))
plot(h, yaxt="n", xaxt="n", cex.main=0.9, col=c("white","black")[ccat], xlim = c(0,maxx), ylim = c(0,maxy), main="papC",xlab="", ylab="")
axis(side=2, at=seq(from=0,to=maxy,by=1000), labels=as.character(seq(from=0,to=maxy,by=1000)), cex.axis = 0.7)
mtext("Frequency", side=2, line=2.2, cex=0.9)
axis(side=1, at=seq(from=0,to=maxx,by=10), labels=as.character(seq(from=0,to=maxx,by=10)), cex.axis = 0.7)
mtext(expression(paste(chi^2)), side=1, line=2.3,cex=1)


xx<-9.619
yy<-0
new.x<-9.619
new.y<-1500
Text<-expression(paste(D^2)==9.619)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=9.619, y=1700, label=Text, cex=0.65)
##############################################################################
##############################################################################
##########################  sat  ################################

v_sat_X<-c(113, 15, 46, 77)
v_sat_pars<-c("L1"=6.003,"L2"=164.8,"L3"=141.0,"L4"= 0,"q12"=0,"q21"=299.9,"q34"=124.1,"q43"=17,"m32"=44.77,"m42"=0)
res<-f_CLTtest(pars=v_sat_pars,v_X=v_sat_X)
Sigma_b<-res$Corr3_16$Sigma_b
v_prop_mean<-res$Corr3_16$v_prop_mean
type_drop<-res$Corr3_16$type_drop
q<-res$Corr3_16$q
res


## order of params L1, L2, L3, L4, mu1, mu2, mu3, mu4, q12, q13(0), q14(0), q21, q23(0), q24(0), q31(0), q32(0), q34, q41(0), q42, q43
v_sat_pars_diversitree<-c(6.003,164.8,141,0,
                          0,0,0,0,
                          0,0,0,
                          299.9,0,0,
                          0,44.77,124.1,
                          0,0,17)
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
sat #6733
subset(sat, sat> 6.732) #741
subset(sat, sat>11.345) #247
max(sat) #32.4355

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


xx<-6.732
yy<-0
new.x<-6.732
new.y<-1500
Text<-expression(paste(D^2)==6.732)
arrows(xx, yy, new.x, new.y, length = 0.07,lwd=1,lty=1,cex=2)
text(x=8, y=1700, label=Text, cex=0.65)
