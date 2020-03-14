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


#Load the data for each VF file (U1 is 2, U0 is 1, K1 is 3, K0 is 4)
#as an example, I load 'fyuA'
data<-as.matrix(read.csv("strains_fyuA.csv",header=TRUE,row.names=1))
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

## or accoring to the observed frequencies similar to what musse by default does (pg. 25 of diversitree manual):

vtip_counts<-table(data.v)
rootpr<-c("1"=0,"2"=0,vtip_counts[3:4])/sum(vtip_counts[3:4])
names(rootpr)<-as.character(1:length(rootpr))
names(rootpr)<-NULL
rootpr

## and then the values sit in rootpr
## but it is NOT enough to write root.p=rootpr
## as then there are some strange warning, it seems you have to manually copy the values ...

defaults <- alist(root = ROOT.GIVEN, root.p = c(0,0,0.601626, 0.398374))

# Modify likelihood function 
lik <- set.defaults(lik, defaults = defaults)
argnames(lik)
p <- starting.point.musse(tr2, 4)
#fit <- find.mle(lik, p)


# Constraints
lik.0 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q42~0,q24~0,
                    q14~0,
                    q23~0,q41~0
)
fit.0 <- find.mle(lik.0, p[argnames(lik.0)])
round(coef(fit.0), 3)


lik.1 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~q21
)
fit.1 <- find.mle(lik.1, p[argnames(lik.1)])
round(coef(fit.1), 3)

lik.2 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~q43
)
fit.2 <- find.mle(lik.2, p[argnames(lik.2)])
round(coef(fit.2), 3)

lik.3 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~q43, q21~q12
)
fit.3 <- find.mle(lik.3, p[argnames(lik.3)])
round(coef(fit.3), 3)

lik.4 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~0
)
fit.4 <- find.mle(lik.4, p[argnames(lik.4)])
round(coef(fit.4), 3)

lik.5 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q21~0
)
fit.5 <- find.mle(lik.5, p[argnames(lik.5)])
round(coef(fit.5), 3)

lik.6 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q34~0
)
fit.6 <- find.mle(lik.6, p[argnames(lik.6)])
round(coef(fit.6), 3)


lik.7 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q43~0
)
fit.7 <- find.mle(lik.7, p[argnames(lik.7)])
round(coef(fit.7), 3)

lik.8 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~q21, q43~0
)
fit.8 <- find.mle(lik.8, p[argnames(lik.8)])
round(coef(fit.8), 3)

lik.9 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,
                   q12~q21,q34~0
)
fit.9 <- find.mle(lik.9, p[argnames(lik.9)])
round(coef(fit.9), 3)

lik.10 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,q12~0,q34~q43
)
fit.10 <- find.mle(lik.10, p[argnames(lik.10)])
round(coef(fit.10), 3)

lik.11 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,q21~0,q34~q43
)
fit.11 <- find.mle(lik.11, p[argnames(lik.11)])
round(coef(fit.11), 3)


lik.12 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0, lambda1~lambda2
)
fit.12 <- find.mle(lik.12, p[argnames(lik.12)])
round(coef(fit.12), 3)


lik.13 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0, lambda3~lambda4
)
fit.13 <- find.mle(lik.13, p[argnames(lik.13)])
round(coef(fit.13), 3)

lik.14 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,lambda1~lambda2,q12~q21
)
fit.14 <- find.mle(lik.14, p[argnames(lik.14)])
round(coef(fit.14), 3)


lik.15 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,lambda3~lambda4,q34~q43
)
fit.15 <- find.mle(lik.15, p[argnames(lik.15)])
round(coef(fit.15), 3)

lik.16 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,
                   q23~0,q41~0,lambda1~lambda2,
                   lambda3~lambda4,q12~q21,q34~q43
)
fit.16 <- find.mle(lik.16, p[argnames(lik.16)])
round(coef(fit.16), 3)

lik.17 <- constrain(lik, 
                   mu1~0,mu2~0,mu3~0,mu4~0,
                   q13~0,q31~0,
                   q42~0,q24~0,
                   q14~0,q41~0,
                   q23~0,
                   lambda1~0
)
fit.17 <- find.mle(lik.17, p[argnames(lik.17)])
round(coef(fit.17), 3)

lik.18 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q42~0,q24~0,
                    q14~0,q41~0,
                    q23~0,
                    lambda2~0
)
fit.18 <- find.mle(lik.18, p[argnames(lik.18)])
round(coef(fit.18), 3)

lik.19 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q42~0,q24~0,
                    q14~0,q41~0,
                    q23~0,
                    lambda3~0
)
fit.19 <- find.mle(lik.19, p[argnames(lik.19)])
round(coef(fit.19), 3)


lik.20 <- constrain(lik, 
                    mu1~0,mu2~0,mu3~0,mu4~0,
                    q13~0,q31~0,
                    q42~0,q24~0,
                    q14~0,q41~0,
                    q23~0,
                    lambda4~0
)
fit.20 <- find.mle(lik.20, p[argnames(lik.20)])
round(coef(fit.20), 3)

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

###############################################################

