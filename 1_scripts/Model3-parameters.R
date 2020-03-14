library(diversitree)
library(ape)


#Load the tree file
tr <- read.tree("tree-CGG+GTG.txt")
plot(tr)
plot.phylo(tr,cex=0.65,no.margin=TRUE)
plot(tr, show.tip.label=FALSE)

#Make this tree ultrametric
tr2 <- chronos(tr)
plot(tr2)
plot(tr2, show.tip.label=FALSE)


#Load the data file (U1 is 2, U0 is 1, K1 is 3, K0 is 4)
#as an example, I load fimG
data<-as.matrix(read.csv("strains_fimG.csv",header=TRUE,row.names=1))
data


data.v<-data[,1]
names(data.v)<-row.names(data)
tr2$tip.state<-data.v
table(data.v)


# Make likelihood function 
lik<-make.musse(tree=tr2,states=tr2$tip.state,k=4)

## the root state setting has to be done due to assumption (F5)

## you may choose how to set the root state probabilties either flat:
# defaults <- alist(root = ROOT.GIVEN, root.p = c(0,0,0.5, 0.5))

## or accoring to the observed frequencies similar to what musse by default does (p. 25 of diversitree manual):

vtip_counts<-table(data.v)
rootpr<-c("1"=0,"2"=0,vtip_counts[3:4])/sum(vtip_counts[3:4])
names(rootpr)<-as.character(1:length(rootpr))
names(rootpr)<-NULL
rootpr

## and then the values sit in rootpr
## but it is NOT enough to write root.p=rootpr
## as then there are some strange warning, it seems you have to manually copy the values ...

defaults <- alist(root = ROOT.GIVEN, root.p = c(0,0,0.97560976, 0.02439024))

# Modify likelihood function 
lik <- set.defaults(lik, defaults = defaults)
argnames(lik)
p <- starting.point.musse(tr2, 4)
#fit <- find.mle(lik, p)


# Constraints
lik.0 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0
)
fit.0 <- find.mle(lik.0, p[argnames(lik.0)])
round(coef(fit.0), 3)


lik.1 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~q21
)
fit.1 <- find.mle(lik.1, p[argnames(lik.1)])
round(coef(fit.1), 3)


lik.2 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~q43
)
fit.2 <- find.mle(lik.2, p[argnames(lik.2)])
round(coef(fit.2), 3)

lik.3 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~q43, q21~q12
)
fit.3 <- find.mle(lik.3, p[argnames(lik.3)])
round(coef(fit.3), 3)


lik.4 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~0
)
fit.4 <- find.mle(lik.4, p[argnames(lik.4)])
round(coef(fit.4), 3)


lik.5 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q21~0
)
fit.5 <- find.mle(lik.5, p[argnames(lik.5)])
round(coef(fit.5), 3)

lik.6 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~0
)
fit.6 <- find.mle(lik.6, p[argnames(lik.6)])
round(coef(fit.6), 3)


lik.7 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q43~0
)
fit.7 <- find.mle(lik.7, p[argnames(lik.7)])
round(coef(fit.7), 3)


lik.8 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~q21, q43~0
)
fit.8 <- find.mle(lik.8, p[argnames(lik.8)])

lik.9 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~q21,q34~0
)
fit.9 <- find.mle(lik.9, p[argnames(lik.9)])
round(coef(fit.9), 3)

lik.10 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,q12~0,q34~q43
)
fit.10 <- find.mle(lik.10, p[argnames(lik.10)])
round(coef(fit.10), 3)


lik.11 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,q21~0,q34~q43
)
fit.11 <- find.mle(lik.11, p[argnames(lik.11)])
round(coef(fit.11), 3)

lik.12 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0, lambda1~lambda2
)
fit.12 <- find.mle(lik.12, p[argnames(lik.12)])
round(coef(fit.12), 3)


lik.13 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0, lambda3~lambda4
)
fit.13 <- find.mle(lik.13, p[argnames(lik.13)])
round(coef(fit.13), 3)


lik.14 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,lambda1~lambda2,q12~q21
)
fit.14 <- find.mle(lik.14, p[argnames(lik.14)])
round(coef(fit.14), 3)


lik.15 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,lambda3~lambda4,q34~q43
)
fit.15 <- find.mle(lik.15, p[argnames(lik.15)])
round(coef(fit.15), 3)

lik.16 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,lambda1~lambda2,
                   lambda3~lambda4,q12~q21,q34~q43
)
fit.16 <- find.mle(lik.16, p[argnames(lik.16)])
round(coef(fit.16), 3)

lik.17 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   lambda1~0
)
fit.17 <- find.mle(lik.17, p[argnames(lik.17)])
round(coef(fit.17), 3)

lik.18 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,
                    lambda2~0
)
fit.18 <- find.mle(lik.18, p[argnames(lik.18)])
round(coef(fit.18), 3)

lik.19 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,
                    lambda3~0
)
fit.19 <- find.mle(lik.19, p[argnames(lik.19)])
round(coef(fit.19), 3)

lik.20 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,
                    lambda4~0
)
fit.20 <- find.mle(lik.20, p[argnames(lik.20)])
round(coef(fit.20), 3)

lik.21 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q42~q32
)
fit.21 <- find.mle(lik.21, p[argnames(lik.21)])
round(coef(fit.21), 3)

lik.22 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~q21,q42~q32
)
fit.22 <- find.mle(lik.22, p[argnames(lik.22)])
round(coef(fit.22), 3)

lik.23 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~q43,q42~q32
)
fit.23 <- find.mle(lik.23, p[argnames(lik.23)])
round(coef(fit.23), 3)

lik.24 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~q43, q21~q12,q42~q32
)
fit.24 <- find.mle(lik.24, p[argnames(lik.24)])
round(coef(fit.24), 3)

lik.25 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~0,q42~q32
)
fit.25 <- find.mle(lik.25, p[argnames(lik.25)])
round(coef(fit.25), 3)

lik.26 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q21~0,q42~q32
)
fit.26 <- find.mle(lik.26, p[argnames(lik.26)])
round(coef(fit.26), 3)

lik.27 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~0,q42~q32
)
fit.27 <- find.mle(lik.27, p[argnames(lik.27)])
round(coef(fit.27), 3)

lik.28 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q43~0,q42~q32
)
fit.28 <- find.mle(lik.28, p[argnames(lik.28)])
round(coef(fit.28), 3)

lik.29 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~q21, q43~0,q42~q32
)
fit.29 <- find.mle(lik.29, p[argnames(lik.29)])
round(coef(fit.29), 3)

lik.30 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,
                    q12~q21,q34~0,q42~q32
)
fit.30 <- find.mle(lik.30, p[argnames(lik.30)])
round(coef(fit.30), 3)

lik.31 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,q12~0,q34~q43,q42~q32
)
fit.31 <- find.mle(lik.31, p[argnames(lik.31)])
round(coef(fit.31), 3)

lik.32 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,q21~0,q34~q43,q42~q32
)
fit.32 <- find.mle(lik.32, p[argnames(lik.32)])
round(coef(fit.32), 3)

lik.33 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0, lambda1~lambda2,q42~q32
)
fit.33 <- find.mle(lik.33, p[argnames(lik.33)])
round(coef(fit.33), 3)

lik.34 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0, lambda3~lambda4,q42~q32
)
fit.34 <- find.mle(lik.34, p[argnames(lik.34)])
round(coef(fit.34), 3)

lik.35 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,lambda1~lambda2,q12~q21,q42~q32
)
fit.35 <- find.mle(lik.35, p[argnames(lik.35)])
round(coef(fit.35), 3)

lik.36 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,lambda3~lambda4,q34~q43,q42~q32
)
fit.36 <- find.mle(lik.36, p[argnames(lik.36)])
round(coef(fit.36), 3)

lik.37 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,lambda1~lambda2,
                    lambda3~lambda4,q12~q21,q34~q43,q42~q32
)
fit.37 <- find.mle(lik.37, p[argnames(lik.37)])

lik.38 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,
                    lambda1~0,q42~q32
)
fit.38 <- find.mle(lik.38, p[argnames(lik.38)])
round(coef(fit.38), 3)

lik.39 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,
                    lambda2~0,q42~q32
)
fit.39 <- find.mle(lik.39, p[argnames(lik.39)])
round(coef(fit.39), 3)

lik.40 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,
                    lambda3~0,q42~q32
)
fit.40 <- find.mle(lik.40, p[argnames(lik.40)])
round(coef(fit.40), 3)

lik.41 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q24~0,
                    q14~0,
                    q23~0,q41~0,
                    lambda4~0,q42~q32
)
fit.41 <- find.mle(lik.41, p[argnames(lik.41)])
round(coef(fit.41), 3)


##BIC
AIC(fit.0,k=log(251))
AIC(fit.1,k=log(251))
AIC(fit.2,k=log(251))
AIC(fit.3,k=log(251))
AIC(fit.4,k=log(251))
AIC(fit.5,k=log(251))
AIC(fit.6,k=log(251))
AIC(fit.7,k=log(251))
AIC(fit.8,k=log(251))
AIC(fit.9,k=log(251))
AIC(fit.10,k=log(251))
AIC(fit.11,k=log(251))
AIC(fit.12,k=log(251))
AIC(fit.13,k=log(251))
AIC(fit.14,k=log(251))
AIC(fit.15,k=log(251))
AIC(fit.16,k=log(251))
AIC(fit.17,k=log(251))
AIC(fit.18,k=log(251))
AIC(fit.19,k=log(251))
AIC(fit.20,k=log(251))
AIC(fit.21,k=log(251))
AIC(fit.22,k=log(251))
AIC(fit.23,k=log(251))
AIC(fit.24,k=log(251))
AIC(fit.25,k=log(251))
AIC(fit.26,k=log(251))
AIC(fit.27,k=log(251))
AIC(fit.28,k=log(251))
AIC(fit.29,k=log(251))
AIC(fit.30,k=log(251))
AIC(fit.31,k=log(251))
AIC(fit.32,k=log(251))
AIC(fit.33,k=log(251))
AIC(fit.34,k=log(251))
AIC(fit.35,k=log(251))
AIC(fit.36,k=log(251))
AIC(fit.37,k=log(251))
AIC(fit.38,k=log(251))
AIC(fit.39,k=log(251))
AIC(fit.40,k=log(251))
AIC(fit.41,k=log(251))
###############################################################
