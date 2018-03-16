### R code from vignette source 'IntroSyn.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Initial R Code
###################################################
options(prompt = "R> ", continue = "+  ", width = 80, useFancyQuotes = FALSE)
options(SweaveHooks=list(fig=function() par(mar=c(4,4,4,1.5)+0.1,cex.axis=1,cex.lab=1,cex.main=1)))
options(SweaveSyntax="SweaveSyntaxNoweb")
set.seed(1234)
library(synbreed)
library(synbreedData)


###################################################
### code chunk number 2: Call to create an object of class gpData (eval = FALSE)
###################################################
## gp <- create.gpData(covar,pheno,geno,map,pedigree)


###################################################
### code chunk number 3: Load the maize data set
###################################################
library(synbreed)
library(synbreedData)
data(maize)


###################################################
### code chunk number 4: IntroSyn.Rnw:387-388
###################################################
summary(maize)


###################################################
### code chunk number 5: maize data: Plots for the genetic map and a histogram of phenotypic values
###################################################
pdf("figs/genMapMaize.pdf")
plotGenMap(maize,ylab="pos[cM]")
dev.off()
pdf("figs/maizeHist.pdf")
hist(maize$pheno[,1,],xlab="yield [dt/ha]",main="",freq=FALSE)
dev.off()


###################################################
### code chunk number 6: Load the mice data set
###################################################
data("mice")


###################################################
### code chunk number 7: mice data: Plot of the genetic map
###################################################
pdf("figs/genMapMice.pdf")
plotGenMap(mice,TRUE,FALSE,ylab="pos [cM]")
dev.off()


###################################################
### code chunk number 8: mice data: Plot scatterplot matrix of the phenotype
###################################################
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="lightgrey")
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y,use="pairwise.complete.obs")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * (0.5+r/2))
}

mypoints <- function(x,y, ...) points(x,y,col=hsv(alpha=0.25,v=0))
pdf("figs/miceHist.pdf")
pairs(mice$pheno,diag.panel = panel.hist,upper.panel=mypoints,lower.panel=panel.cor)
dev.off()


###################################################
### code chunk number 9: Load the cattle data set
###################################################
data("cattle")


###################################################
### code chunk number 10: mice data: Plot of the genetic map
###################################################
pdf("figs/genMapCattle.pdf")
plotGenMap(cattle,TRUE,FALSE,ylab="pos [Mb]")
dev.off()


###################################################
### code chunk number 11: cattle data: Plot scatterplot matrix of the phenotype
###################################################
pdf("figs/cattleHist.pdf")
pairs(cattle$pheno,diag.panel = panel.hist,upper.panel=mypoints,lower.panel=panel.cor)
dev.off()


###################################################
### code chunk number 12: Recode maize data
###################################################
maizeC <- codeGeno(maize)


###################################################
### code chunk number 13: Recode mice data (eval = FALSE)
###################################################
## # simple recoding of alleles
## miceC <- codeGeno(mice,label.heter="alleleCoding")


###################################################
### code chunk number 14: Recode and impute mice data (eval = FALSE)
###################################################
## # recoding of mice data and imputing of missing values by family structure
## # discarding markers with maf < 0.05 and nmiss > 0.01
## miceC <- codeGeno(mice,label.heter=is.heter,impute=TRUE,impute.type="random",maf=0.05,nmiss=0.01,verbose=TRUE)


###################################################
### code chunk number 15: Pairwise LD for the maize data (eval = FALSE)
###################################################
## maizeLD <- pairwiseLD(maizeC,chr=1,type="data.frame")


###################################################
### code chunk number 16: LD decay scatterplot (eval = FALSE)
###################################################
## LDDist(maizeLD,type="p",xlab="dist [cM]",pch=19,col=hsv(alpha=0.075,v=0))


###################################################
### code chunk number 17: LD decay stacked histogram (eval = FALSE)
###################################################
## LDDist(maizeLD,type="bars",breaks=list(dist=c(0,25,50,75,200),
## r2=c(1,0.5,0.3,0.2,0.1,0.05,0)),xlab="dist [cM]")


###################################################
### code chunk number 18: Additive numerator relationship matrix for the maize data (eval = FALSE)
###################################################
## A <- kin(maize,ret="add",DH=maize$covar$DH)


###################################################
### code chunk number 19: Summary method for class relationshipMatrix (eval = FALSE)
###################################################
## summary(A)


###################################################
### code chunk number 20: Realized relationship matrix for the maize data (eval = FALSE)
###################################################
## U <- kin(maizeC,ret="realized")
## summary(U)


###################################################
### code chunk number 21: Heatmap visualization of the expected and realized relationship matrices (eval = FALSE)
###################################################
## plot(A[maize$covar$genotyped,maize$covar$genotyped])
## plot(U)


###################################################
### code chunk number 22: Realized relationship matrix for the mice data (eval = FALSE)
###################################################
## UM <- kin(miceC,ret="realized")


###################################################
### code chunk number 23: Genomic BLUP for the mice data using trait weight (eval = FALSE)
###################################################
## miceGBLUP <- gpMod(miceC,model="BLUP",kin=UM,trait="weight")


###################################################
### code chunk number 24: Summary of the model (eval = FALSE)
###################################################
## summary(miceGBLUP)


###################################################
### code chunk number 25: Genomic BLUP for the mice data using trait weight (eval = FALSE)
###################################################
## miceRRBLUP <- gpMod(miceC,model="BLUP",kin=UM,trait="weight",markerEffects=TRUE)


###################################################
### code chunk number 26: IntroSyn.Rnw:675-682
###################################################
# 40 = 8 * (3+2)
# h2 = 0.74 (Valdar 2006)
# 52.3  = sqrt(2*(1-h2)/h2*sum(colMeans(miceC$geno)^2))
# lambda = 1:100
# plot(dgamma(lambda^2,shape=0.8,rate=1e-4)*lambda*2)
# abline(v=52.3)
prior <- list(varE=list(df=3,S=40),lambda = list(shape=0.8,rate=1e-4,value=52,type='random'))


###################################################
### code chunk number 27: Bayesian Lasso model for the mice data (eval = FALSE)
###################################################
## miceModBL  <- gpMod(miceC,model="BL",trait="weight",
## prior=prior,nIter=12000,burnIn=2000,thin=10)


###################################################
### code chunk number 28: Summary of the model (eval = FALSE)
###################################################
## summary(miceModBL)


###################################################
### code chunk number 29: Manhattanplot for the marker effects in the mice data (eval = FALSE)
###################################################
## pdf("figs/manhattanPlot.pdf",width=11,height=5)
## manhattanPlot(abs(miceModBL$m),miceC,ylab="|marker effect|")
## dev.off()


###################################################
### code chunk number 30: Set up a gpData for the prediction set (eval = FALSE)
###################################################
## unphenotyped <- dimnames(mice$pheno)[[1]][is.na(mice$pheno[,1,])]
## phenotyped <- mice$covar$id[!mice$covar$id %in% unphenotyped]
## predSet <- discard.individuals(miceC,phenotyped)


###################################################
### code chunk number 31: Predict their genetic performance (eval = FALSE)
###################################################
## predict(miceModBL,predSet)


###################################################
### code chunk number 32: Cross-validation of the GBLUP model with the mice data (eval = FALSE)
###################################################
## cv.mice <- crossVal(miceC,cov.matrix=list(UM),k=2,Rep=10,Seed=123,
##              sampling="random",varComp=miceGBLUP$fit$sigma,VC.est="commit")


###################################################
### code chunk number 33: Load previous results to save time
###################################################
load("cv.mice.Rdata")
cv.mice <- cv.mice2


###################################################
### code chunk number 34: Summary for the cross-validation
###################################################
summary(cv.mice)


###################################################
### code chunk number 35: Models for the maize data (eval = FALSE)
###################################################
## # animal model
## PBLUP <- gpMod(maizeC,model="BLUP",kin=A/2)
## # G-BLUP
## GBLUP <- gpMod(maizeC,model="BLUP",kin=U/2)
## # Bayesian Lass
## prior <- list(varE=list(df=3,S=35),lambda = list(shape=0.52,rate=1e-4,value=20,type='random'))
## modBL <- gpMod(maizeC,model="BL",prior=prior,nIter=6000,burnIn=1000,thin=5)


###################################################
### code chunk number 36: Extract true breeding values from the maize data object
###################################################
tbv <- maize$covar$tbv[maize$covar$genotyped]


