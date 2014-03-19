#############################################
## make_maps.R
#
# purpose:	compute genetic maps, variously subdivided, from CC G2:F1 population
# author:		GAC 1 May 2012
#	modified:	APM 8 Oct 2013
#
#############################################

## R libraries
library(ggplot2)
library(plyr)

## load functions library
source("./helper_functions.R") 
ls()

## define colors for CC founders
CC.colors <- c(AJ="#F0F000",NOD="#1010F0",B6="#404040",S129="#F08080",
  WSB="#9000E0",PWK="#F00000",NZO="#00A0F0",CAST="#00A000")
CC.colors.alt <- c(AJ="#DAA520",NOD="#1010F0",B6="#404040",S129="#F08080",
  WSB="#9000E0",PWK="#F00000",NZO="#00A0F0",CAST="#00A000")  
color.sex <- c(Male="blue",Female="red")

## names for factor levels
StrainNames <- c("S129","AJ","B6","CAST","NOD","NZO","PWK","WSB")
SubSpecNames <- c("CAST","DOM","MUS","UNK")
SubSpecPairNames <- c("CAST.DOM","CAST.MUS","DOM.DOM","DOM.MUS","UNK")
Prdm9Alleles <- c("A","B","C","D")
# A is the CAST allele w/ 11 fingers
# B is B6/S129/AJ/NZO  w/ 12 fingers
# C is WSB/NOD  w/13 fingers
# D is the PWK allele with 14 fingers
Prdm9Names <- c("A","AB","AC","AD","B","BC","BD","C","CD","D") # expanded to include ambiguity codes
MeiosisNames <- c("MGM","MGP","PGM","PGP","M","P") # MGM, PGM and M are the female meioses
MeiosisNames.new <- c("MGM","MGP","PGM","PGP","Mf","Pf","Mm","Pm") # extended version distinguishes meioses that give rise to each G2F1 animal
GenNames <- c("G1", "G2") # generation of the meiosis event
SexNames=c("Male","Female")


#############################################
## load and annotate G2F1 recombination events data
# NB: this is the same data as in FileS1.csv, but with different column headings

Events.all <- read.csv("../data/AllEvents6.csv", header = TRUE, stringsAsFactors = TRUE, row.names = NULL)
dim(Events.all)
# [1] 25038    27

# remove one copy of each shared event from Events.all dataframe
Events <- Events.all[!(Events.all[,"Shared"]=="Shared" & Events.all[,"Sex"]=="f"),]
dim(Events)
# [1] 21993    27

# in this script we are only counting shared events once and will not need Events.all
rm(Events.all)
ls()

# break out the G2 meioses according to which G2F1 animal they give rise to
Type.new <- ifelse(Events$Type=="M" & Events$Sex=="f", "Mf", Events$Type)
Type.new <- ifelse(Events$Type=="M" & Events$Sex=="m", "Mm", Type.new)
Type.new <- ifelse(Events$Type=="P" & Events$Sex=="f", "Pf", Type.new)
Type.new <- ifelse(Events$Type=="P" & Events$Sex=="m", "Pm", Type.new)
Type.new <- ifelse(Events$Type=="PGM", "PGM", Type.new)
Type.new <- ifelse(Events$Type=="MGM", "MGM", Type.new)
Type.new <- ifelse(Events$Type=="PGP", "PGP", Type.new)
Type.new <- ifelse(Events$Type=="MGP", "MGP", Type.new)
Events <- cbind(Events, Type.new)
# there are eight different meioses event in each funnel

# assign generation
Gen <- rep(NA, nrow(Events))
Gen[ grep("G", as.character(Events$Type.new)) ] <- "G1"
Gen[ is.na(Gen) ] <- "G2"
Events <- cbind(Events, Gen)

## count funnels and create indices to select unique records corresponding to meioses or funnels
# how many unique funnels?
(n.funnels <- length(unique(Events[,"Funnel_ID"])))
# there are 237 unique funnels

## compute an index to capture each meiosis once
indx.m <- which(!duplicated(Events[,c("Funnel_ID","Type.new")]))
length(indx.m)
# there are 1896 = 237 * 8 distinct meioses

## compute an index to capture each funnel once
indx.f <- which(!duplicated(Events[,"Funnel_ID"]))
length(indx.f)

## check that there are 8 meioses per funnel
length(indx.m)/length(indx.f)

## add a variable to track the sex of the meiosis events
Events <- transform(Events, 
	Sex.M = ifelse(Type=="M"|Type=="MGM"|Type=="PGM","Female","Male"))

## number of events per sex (autosome)
table(subset(Events,Chr<20)$Sex.M)
#Female   Male 
# 11241  10127 

## check that number of meioses is 4 per sex
table(Events[indx.m,]$Sex.M)/4
# 237

## compute indicator variable for strains in each event junction
Events <- transform(Events,
 S129 = (Dist_Founder=="S129" | Prox_Founder=="S129"),
 AJ = (Dist_Founder=="AJ" | Prox_Founder=="AJ"),
 B6 = (Dist_Founder=="B6" | Prox_Founder=="B6"),
 CAST = (Dist_Founder=="CAST" | Prox_Founder=="CAST"),
 NOD = (Dist_Founder=="NOD" | Prox_Founder=="NOD"), 
 NZO = (Dist_Founder=="NZO" | Prox_Founder=="NZO"),
 PWK = (Dist_Founder=="PWK" | Prox_Founder=="PWK"),
 WSB = (Dist_Founder=="WSB" | Prox_Founder=="WSB")
)
# the indicator variables are handy because each event involves two strains
# indicator is TRUE if either the strain is either distal or proximal  

## how many junctions does each strain contribute to?
apply(Events[,StrainNames],2,sum)
#S129   AJ   B6 CAST  NOD  NZO  PWK  WSB 
#5573 5434 5555 5470 5572 5537 5379 5466 

## check that the sum is equal to total number of events
sum(apply(Events[,StrainNames],2,sum))/2
# 21993

## tally all of the junctions by strain (ordered)
table(Events$Dist_Founder, Events$Prox_Founder)
#        AJ  B6 CAST NOD NZO PWK S129 WSB
#  AJ     0 407  280 276 503 435  371 417
#  B6   409   0  390 370 513 392  289 398
#  CAST 277 407    0 486 331 470  422 348
#  NOD  306 326  515   0 416 429  407 418
#  NZO  507 516  360 374   0 317  361 318
#  PWK  415 408  431 443 332   0  346 306
#  S129 420 303  446 390 358 373    0 535
#  WSB  411 427  307 416 331 282  552   0

## compute the subspecific type of each junction
# there are more than four possible combinations
# but we will focus on four primary types and call the rest unknown
levels(Events$Prox_Origin)
levels(Events$Dist_Origin)

Events <- transform(Events,
	Junction.SubType = 
		ifelse((Prox_Origin=="Castaneous"&Dist_Origin=="Domesticus")|
			(Dist_Origin=="Castaneous"&Prox_Origin=="Domesticus"),"CAST.DOM",
		ifelse((Prox_Origin=="Castaneous"&Dist_Origin=="Musculus")|
			(Dist_Origin=="Castaneous"&Prox_Origin=="Musculus"),"CAST.MUS",
		ifelse((Prox_Origin=="Domesticus"&Dist_Origin=="Domesticus"),"DOM.DOM",
		ifelse((Prox_Origin=="Musculus"&Dist_Origin=="Domesticus")|
			(Dist_Origin=="Musculus"&Prox_Origin=="Domesticus"),"DOM.MUS","UNK")))) )

## tabulate junction subspecies types
table(Events$Junction.SubType)

## check that it matches total number of events
sum(table(Events$Junction.SubType))
#  21987 !!! 6 events dropped due to incomplete data

## tally all junctions by subspecies, ordered
table(Events$Prox_Origin, Events$Dist_Origin)
#             Castaneous Domesticus Multiple Musculus Unknown
#  Castaneous         19       1852        2      476      18
#  Domesticus       1828      11886       15     2267      94
#  Multiple            5         24      117        4       0
#  Musculus          507       2249        5      457      21
#  Unknown             5        104        0       32       0

## check total number of events
tmp <- table(Events$Prox_Origin, Events$Dist_Origin)
sum(apply(tmp,1,"sum"))
#  21987  this sum is also short 6 events

## assign PRDM9 allele types at each event
# convert the quarter funnel encoding to a Prdm9 genotype at each meiosis
Events <- transform(Events,
PRDM9_Q1A = factor(levels(PRDM9_Q1A)[PRDM9_Q1A], levels=Prdm9Names), 
PRDM9_Q1B = factor(levels(PRDM9_Q1B)[PRDM9_Q1B], levels=Prdm9Names), 
PRDM9_Q2A = factor(levels(PRDM9_Q2A)[PRDM9_Q2A], levels=Prdm9Names), 
PRDM9_Q2B = factor(levels(PRDM9_Q2B)[PRDM9_Q2B], levels=Prdm9Names), 
PRDM9_Q3A = factor(levels(PRDM9_Q3A)[PRDM9_Q3A], levels=Prdm9Names), 
PRDM9_Q3B = factor(levels(PRDM9_Q3B)[PRDM9_Q3B], levels=Prdm9Names), 
PRDM9_Q4A = factor(levels(PRDM9_Q4A)[PRDM9_Q4A], levels=Prdm9Names), 
PRDM9_Q4B = factor(levels(PRDM9_Q4B)[PRDM9_Q4B], levels=Prdm9Names)) 
for(i in 16:27){
	print(table(Events[,i]))
}
# note that double letters indicate ambiguity in assignment, NOT genotypes
# eg, AB = A or B allele

## PRDM9L is the maternal allele and PRDM9R is the paternal allele
Events <- transform(Events, 
	Prdm9L = as.factor(ifelse(Type=="MGM",PRDM9_Q1A,
		ifelse(Type=="MGP",PRDM9_Q2A,
			ifelse(Type=="PGM",PRDM9_Q3A,
				ifelse(Type=="PGP",PRDM9_Q4A,
					ifelse(Type=="M",PRDM9_F1LA,PRDM9_F1RA)))))),
	Prdm9R = as.factor(ifelse(Type=="MGM",PRDM9_Q1B,
		ifelse(Type=="MGP",PRDM9_Q2B,
			ifelse(Type=="PGM",PRDM9_Q3B,
				ifelse(Type=="PGP",PRDM9_Q4B,
					ifelse(Type=="M",PRDM9_F1LB,PRDM9_F1RB)))))) )

levels(Events$Prdm9L) <- Prdm9Names
levels(Events$Prdm9R) <- Prdm9Names
table(Events$Prdm9L, Events$Prdm9R)

## how many meioses are associated with each Prdm9 genotype
table(Events[indx.m,]$Prdm9L, Events[indx.m,]$Prdm9R)[c(1,5,8,10),c(1,5,8,10)]
#      A   B   C   D
#  A   0 114  68  41
#  B 109 331 236 140
#  C  53 208  44  52
#  D  24  77  57   0

## how many events are associated with each Prdm9 genotype (ignoring ambiguous cases)
table(Events$Prdm9L, Events$Prdm9R)[c(1,5,8,10),c(1,5,8,10)]
#       A    B    C    D
#  A    0 1355  757  469
#  B 1215 3628 2586 1662
#  C  641 2276  481  600
#  D  283  812  635    0

# breakdown by generation
table(Events$Prdm9L, Events$Prdm9R, Events$Gen)[c(1,5,8,10),c(1,5,8,10),]

## create subsepecific Prdm9 types
# CAS is A,  DOM is B or C,  MUS is D, UNK is ambiguous subspecies
Events <- transform(Events,
	Prdm9L.Sub = ifelse(Prdm9L=="A","CAST",
					ifelse(Prdm9L=="B"|Prdm9L=="BC"|Prdm9L=="C","DOM",
						ifelse(Prdm9L=="D","MUS","UNK"))),
	Prdm9R.Sub = ifelse(Prdm9R=="A","CAST",
					ifelse(Prdm9R=="B"|Prdm9R=="BC"|Prdm9R=="C","DOM",
						ifelse(Prdm9R=="D","MUS","UNK")))	)					

levels(Events$Prdm9L.Sub) <- c("CAST","DOM","MUS","UNK")
levels(Events$Prdm9R.Sub) <- c("CAST","DOM","MUS","UNK")						

## tabulate Prdm9 subspecific genotypes for all events
table(Events$Prdm9L.Sub, Events$Prdm9R.Sub)[1:3,1:3]
#       CAST   DOM  MUS
#  CAST    0  2169  469
#  DOM  2032 10416 2311
#  MUS   283  1682    0

## how many meioses are associated with each Prdm9 subspecific origin
# this table is our demoninator for estimating maps by Prdm9 subspecific type
table(Events[indx.m,]$Prdm9L.Sub, Events[indx.m,]$Prdm9R.Sub)[1:3, 1:3]

# breakdown by generation
table(Events$Prdm9L.Sub, Events$Prdm9R.Sub, Events$Gen)[1:3,1:3,]

# breakdown by Sex.M (sex of the parent in which the meiosis occurred)
table(Events$Prdm9L.Sub, Events$Prdm9R.Sub, Events$Sex.M)[1:3,1:3,]

# breakdown by Gen & Sex.M
table(Events$Prdm9L.Sub, Events$Prdm9R.Sub, Events$Gen, Events$Sex.M)[1:3,1:3,,]

table(Events$Gen, Events$Prdm9L.Sub)
table(Events$Gen, Events$Prdm9R.Sub)



#############################################
## Summary stats and analysis of G2F1 events

## compute parameters of each chromosome
# position of first and last events, length (Mb) of chromosome spanned by events, 
# total number of events, density of events per 100Mb per meiosis
# note correction factor for X chr is 3.75
ChrLen <- cbind(ddply(Events[,c("Start","Chr")],.(Chr),numcolwise(min)),
	ddply(Events[,c("End","Chr")],.(Chr),numcolwise(max)))[,-3]
ChrLen <- transform(ChrLen, 	
	n.Events = ddply(Events, .(Chr), dim)[,2],
	Size.Mb=(End-Start)/1000000)	
ChrLen <- transform(ChrLen, Size.cM = ifelse(Chr<20, 
	100*n.Events/(n.funnels*7),100*n.Events/(n.funnels*3.75)))
ChrLen <- transform(ChrLen,cMperMb = Size.cM/Size.Mb, MbpercM = Size.Mb/Size.cM)
round(ChrLen,2)

## plot physical vs genetic size of chromosomes
qplot(Size.Mb, Size.cM, color=Chr, data=ChrLen, main="G2F1 Chr Size") + 
	xlab("Physical Size (Mb)") + ylab("Genetic Size (cM)") +
	# scale_color_gradientn(colour=rainbow(20), breaks=1:20) +
	geom_smooth(se=F) + geom_point(size=3) +theme_bw()

## helper function to compute scaling factors for cumulative maps
get.scaling.factors <- function(Events, StrainNames = StrainNames, MeiosisNames = MeiosisNames, ...) {
	## count occurence of each strain by each funnel position 8x8 matrix
	Strain.Pos <- NULL
	for(i in 1:8){
		Strain.Pos <- cbind(Strain.Pos, 
												table(substring(unique(Events$Funnel_Order),i,i)))	
	}
	# assign position names A-H
	colnames(Strain.Pos) <- rownames(Strain.Pos)
	# reorder rows and assign strain names
	Strain.Pos <- Strain.Pos[c(3,1,2,6,4,5,7,8),]
	rownames(Strain.Pos) <- StrainNames
	# (show the result)
	# Strain.Pos
	
	# count strain contributions to each meiosis
	Strain.Meiosis <- cbind(Strain.Pos[,"A"]+Strain.Pos[,"B"],Strain.Pos[,"C"]+Strain.Pos[,"D"],
													Strain.Pos[,"E"]+Strain.Pos[,"F"],Strain.Pos[,"G"]+Strain.Pos[,"H"], 
													Strain.Pos[,"A"]+Strain.Pos[,"B"]+Strain.Pos[,"C"]+Strain.Pos[,"D"],
													Strain.Pos[,"E"]+Strain.Pos[,"F"]+Strain.Pos[,"G"]+Strain.Pos[,"H"])
	colnames(Strain.Meiosis) <- MeiosisNames
	# (show the result)
	# Strain.Meiosis
	
	## number of observable meioses per funnel position (see Supplementary Text)
	# first for autosomes...
	Pos.Meiosis.A <- matrix(c(
		3/4,   0,   0,   0, 1, 0,
		3/4,   0,   0,   0, 1, 0,
		0, 3/4,   0,   0, 1, 0,	
		0, 3/4,   0,   0, 1, 0,
		0,   0, 3/4,   0, 0, 1,	
		0,   0, 3/4,   0, 0, 1,	
		0,   0,   0, 3/4, 0, 1,	
		0,   0,   0, 3/4, 0, 1), nrow=8,byrow=T)
	rownames(Pos.Meiosis.A) <- colnames(Strain.Pos)
	colnames(Pos.Meiosis.A) <- MeiosisNames
	# (show the result)
	# Pos.Meiosis.A
	
	# ...and now the X chromosome
	Pos.Meiosis.X <- matrix(c(
		3/4, 0, 0, 0, 1, 0,
		3/4, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 1, 0,	
		0, 0, 0, 0, 1, 0,
		0, 0, 1, 0, 0, 0,	
		0, 0, 1, 0, 0, 0,	
		0, 0, 0, 0, 0, 0,	
		0, 0, 0, 0, 0, 0), nrow=8,byrow=T)
	rownames(Pos.Meiosis.X) <- colnames(Strain.Pos)
	colnames(Pos.Meiosis.X) <- MeiosisNames
	# (show the result)
	# Pos.Meiosis.X
	
	## summation vector for strains
	# note use of 1/2 because two strains in each junction
	Strain.indx <- matrix(rep(1/2,8),nrow=1)
	names(Strain.indx) <- StrainNames
	
	## summation vector for all meioses
	Meiosis.indx <- matrix(rep(1,6),ncol=1)
	names(Meiosis.indx) <- MeiosisNames
	
	return( list(Strain.Pos = Strain.Pos, Strain.Meiosis = Strain.Meiosis,
							 Pos.Meiosis.A = Pos.Meiosis.A, Pos.Meiosis.X = Pos.Meiosis.X,
							 Strain.indx = Strain.indx, Meiosis.indx = Meiosis.indx) )
	
}

my.scaling <- get.scaling.factors(Events, StrainNames, MeiosisNames)
Strain.indx <- my.scaling$Strain.indx
Strain.Pos <- my.scaling$Strain.Pos
Meiosis.indx <- my.scaling$Meiosis.indx
Pos.Meiosis.A <- my.scaling$Pos.Meiosis.A
Pos.Meiosis.X <- my.scaling$Pos.Meiosis.X

## compute some counts as a sanity check
# number of observable events per meiosis for autosomes should be 7
Strain.indx %*% Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx / n.funnels
# strain marginal
Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx / (2 * n.funnels)
# meiosis marginal
Strain.indx %*% Strain.Pos %*% Pos.Meiosis.A  / n.funnels

# number of observable events per meiosis for chrX should be 3.75
Strain.indx %*% Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx / n.funnels
# strain marginal
# NB: strains make different contributions to X due to unequal funnels assignments
Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx / (2 * n.funnels)
# meiosis marginal
Strain.indx %*% Strain.Pos %*% Pos.Meiosis.X  / n.funnels



#############################################
## Summary stats and analysis of map constructed from Heterogeneous Stock, for comparison
# reference: Cox et al (2009) Genetics 182(4): 1335-1344.

## load data from Jax Mouse Map Converter tool
CoxMap <- read.csv(
        "http://cgd.jax.org/mousemapconverter/Revised_HSmap_SNPs.csv",
        colClasses = c(rep("character", 1), rep("numeric", 5)),
        header = TRUE,
        stringsAsFactors = TRUE,
        row.names = NULL)
names(CoxMap) <- c("snpID","Chr","Pos","female_cM","male_cM","average_cM")
# (peek to check data integrity)
head(CoxMap)

## remove the zero markers
CoxMap <- subset(CoxMap, Pos>0)

## shift the maps to place first markers at zero cM
Cox.Start.cM <- ddply(CoxMap,.(Chr),numcolwise(min))[,3:5]
CoxMap$female_cM <- CoxMap$female_cM - Cox.Start.cM[CoxMap$Chr,1]
CoxMap$male_cM <- CoxMap$male_cM - Cox.Start.cM[CoxMap$Chr,2]
CoxMap$average_cM <- CoxMap$average_cM - Cox.Start.cM[CoxMap$Chr,3]

## copy the female_cM into average_cM on X chromosome
CoxMap[CoxMap[,"Chr"]==20,"average_cM"] <- CoxMap[CoxMap[,"Chr"]==20,"female_cM"]

## convert to cumulative map format
CoxMap.avg <- CoxMap[,c("Chr","Pos","average_cM")]
events <- CoxMap.avg$average_cM * 35.46
events[CoxMap.avg$Chr==20] <- events[CoxMap.avg$Chr==20]/2
CoxMap.avg <- cbind(CoxMap.avg, events)
names(CoxMap.avg) <- c("Chr", "Pos", "cM", "cumEvents")

## get sex-specific maps
CoxMap.female <- CoxMap[,c("Chr","Pos","female_cM")]
events <- CoxMap.female$female_cM * 17.73
CoxMap.female <- cbind(CoxMap.female, events)
names(CoxMap.female) <- c("Chr", "Pos", "cM", "cumEvents")
CoxMap.male <- CoxMap[CoxMap[,"Chr"]!=20,c("Chr","Pos","male_cM")]
events <- CoxMap.male$male_cM * 17.73
CoxMap.male <- cbind(CoxMap.male, events)
names(CoxMap.male) <- c("Chr", "Pos", "cM", "cumEvents")

## combine male and female maps
CoxMap.Sex <- rbind(transform(CoxMap.male,Sex="Male"),
                      transform(CoxMap.female,Sex="Female"))
                      
# (peek to check result)
head(CoxMap.avg)
head(CoxMap.Sex)

## (clean up unneeded intermediates)
rm(CoxMap, CoxMap.male, CoxMap.female)
ls()

## compute chromosome sizes in bp, events and cM for Cox maps
ChrLen.Cox <- cbind(ddply(subset(CoxMap.avg, Pos>0),.(Chr),numcolwise(min)),
	ddply(subset(CoxMap.avg, Pos>0),.(Chr),numcolwise(max)))[,c(1,2,6,3,7)]
names(ChrLen.Cox) <- c("Chr","Start.Mb","End.Mb","Start.cM","End.cM")
ChrLen.Cox <- transform(ChrLen.Cox, 
		Size.Mb=(End.Mb-Start.Mb)/1000000, Size.cM=End.cM-Start.cM)
ChrLen.Cox <- transform(ChrLen.Cox, 	
		cMperMb = Size.cM/Size.Mb, MbpercM = Size.Mb/Size.cM)
## remove redundant columns
ChrLen.Cox <- ChrLen.Cox[,-c(4,5)]
# (show result)
round(ChrLen.Cox, 2)

## plot physical vs genetic size of chromosomes
qplot(Size.Mb, Size.cM, color=Chr, data=ChrLen.Cox, main="CoxMap Chr Size") + 
	xlab("Physical Size (Mb)") + ylab("Genetic Size (cM)") +
	scale_color_gradientn(colour=rainbow(20), breaks=1:20) +
	geom_smooth(se=F) + geom_point(size=3) +theme_bw()
	
## compare genetic map lengths of Cox vs G2F1
tmp <- cbind(ChrLen[,c("Chr","Size.cM","MbpercM")], ChrLen.Cox[,c("Size.cM","MbpercM")])
names(tmp) <- c("Chr","Size.cM.G2F1","MbpercM.G2F1","Size.cM.Cox","MbpercM.Cox")

qplot(Size.cM.G2F1, Size.cM.Cox, color=Chr, data=tmp, main="Chr Size" )  + 
	xlab("G2F1 (cM)") + ylab("Cox (cM)") +
	scale_color_gradientn(colour=rainbow(20), breaks=1:20) +
	geom_abline(b=1, a=0) +
	geom_point(size=3) + theme_bw()



###########################################################
##  Compute sex-averaged G2:F1 map

## convert Events to unscaled cumulative events
G2Map.avg <- padMap(ddply(Events, .(Chr), integrateCrossovers, .progress = "text"))
# (same, but on bootstrap samples)
events.bs <- lapply(1:1000, function(i) {
	Events[ sort(sample.int(nrow(Events), replace = TRUE)), ]
})
G2Map.bs <- lapply(events.bs, function(df) {
	padMap(ddply(df, .(Chr), integrateCrossovers))
})

## check that events per chromosome is correct
cbind(ChrLen[,c("Chr","n.Events")],ddply(G2Map.avg,.(Chr),numcolwise(max)))[,c(2,5)]

## scale cumulative map to cM
# count the number of observe meioses in 237 unique funnels x 100 to obtain cM
scalefac.A <- 100 /(Strain.indx %*% Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx)
scalefac.X <- 100 /(Strain.indx %*% Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx)
G2Map.avg <- transform(G2Map.avg, cM = ifelse(Chr<20, 
	cumEvents*scalefac.A, cumEvents*scalefac.X) )
# (do it on bootstrap samples)
G2Map.bs <- lapply(1:length(events.bs), function(i) {
	sf <- get.scaling.factors(events.bs[[1]], StrainNames, MeiosisNames)
	scalefac.A <- 100 /(sf$Strain.indx %*% sf$Strain.Pos %*% sf$Pos.Meiosis.A %*% sf$Meiosis.indx)
	scalefac.X <- 100 /(sf$Strain.indx %*% sf$Strain.Pos %*% sf$Pos.Meiosis.X %*% sf$Meiosis.indx)
	transform(G2Map.bs[[i]], cM = ifelse(Chr<20, cumEvents*scalefac.A, cumEvents*scalefac.X) )
})
# check chr lengths
cbind(ChrLen$Size.cM, ddply(G2Map.avg,.(Chr),numcolwise(max))$cM)

## plot cumulative maps
qplot(Pos, cM, data=G2Map.avg, geom="line") + 
	facet_wrap(facets=~Chr, nrow=5)
qplot(Pos, cM, data=G2Map.avg, geom="line", group=as.factor(Chr), color=Chr) +
	scale_color_gradientn(colours=rainbow(20), breaks=1:20) 



#####################################
## Compute cumulative G2:F1 map by sex

G2Map.Sex <- rbind(
    transform(padMap(ddply(
      Events[apply(outer(Events[,"Type"],c("MGP","PGP","P"),"=="),1,any),],
      	.(Chr),integrateCrossovers)),Sex="Male"),
    transform(padMap(ddply(
      Events[apply(outer(Events[,"Type"],c("MGM","PGM","M"),"=="),1,any),],
      	.(Chr),integrateCrossovers)),Sex="Female"))
      	
## scale cumulative maps to cM units
Meiosis.indx.sex <- matrix(c(0,1,0,1,0,1,1,0,1,0,1,0),ncol=2)
rownames(Meiosis.indx.sex) <- MeiosisNames
colnames(Meiosis.indx.sex) <- c("Male","Female")

scalefac.A <- 100 /Strain.indx %*% Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx.sex
scalefac.X <- 100 /Strain.indx %*% Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx.sex

G2Map.Sex <- transform(G2Map.Sex, 
	cM = ifelse(Chr<20, cumEvents*scalefac.A[Sex],cumEvents*scalefac.X[2] ))
head(G2Map.Sex)

## check that the average chr length across strains is the same as avg chr length
#  chrX will be a little off due to unequal weights
ddply(ddply(G2Map.Sex,.(Chr, Sex),numcolwise(max))[,c(1,5)],.(Chr),numcolwise(mean))[,2] /
	(ddply(G2Map.avg,.(Chr),numcolwise(max))[,4])

## check sex ratio
subset(ddply(G2Map.Sex,.(Chr,Sex),numcolwise(max)),Sex=="Female")[1:19,"cM"] /
	subset(ddply(G2Map.Sex,.(Chr,Sex),numcolwise(max)),Sex=="Male")[,"cM"]

## update chromosome lengths
ChrLen <- transform(ChrLen,
	male_cM = c(ddply(subset(G2Map.Sex,Sex=="Male"),.(Chr, Sex),	
		numcolwise(max))[,"cM"],NA),
	female_cM = ddply(subset(G2Map.Sex,Sex=="Female"),.(Chr, Sex),
		numcolwise(max))[,"cM"])
ChrLen <- transform(ChrLen, Sex_Ratio = female_cM/male_cM)
ChrLen 

## total autosomal map length by Sex 		
ddply(ddply(subset(G2Map.Sex,Chr<20),.(Chr, Sex),numcolwise(max))[,c(2,5)],
	.(Sex),numcolwise(sum))

## plot cumulative maps by sex
qplot(Pos, cM, data=G2Map.Sex, geom="line",color=Sex) + 
	aes(color=Sex) + scale_colour_manual(values=color.sex) +
	facet_wrap(facets=~Chr, nrow=5) 
	 
## plot chr lengths in cM by sex
qplot(male_cM, female_cM, data=subset(ChrLen,Chr<20), colour=Chr) + 
	geom_abline(slope=1,intercept=0,lty=2) +  geom_point(size=3) +
	scale_color_gradientn(colour=rainbow(20), breaks=1:19) + theme_bw()

## plot cumulative maps by sex, 1 page per chromosome, and save it to file
pdf("CumulativeBySex.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
	print(qplot(Pos, cM, data=subset(G2Map.Sex,Chr==i), 
		geom="line", color=Sex,
		main=paste("Chromosome", as.character(i), sep=" ")) +
		aes(color=Sex) + scale_colour_manual(value=color.sex) )
	}
dev.off()



################################################################
## compute density plots (cM/bp) for sex-averaged maps

## compare total events G2F1 vs Cox
total.events <- cbind(ddply(G2Map.avg, .(Chr), numcolwise(max)), 
      ddply(CoxMap.avg, .(Chr), numcolwise(max))[,-1])[,c(1,3,7)]
names(total.events) <- c("Chr","G2F1.Events","Cox.Events")
total.events
	
## prepare to compute density in sliding window: scale = window size; step = offset
# set scale and step size for density computation in kb, 1000=1Mb
dscale = 50
dstep = dscale/5

## G2F1 density
G2Den.avg <- ddply(G2Map.avg, .(Chr), make.density, dstep, dscale)
head(G2Den.avg)
# density in windows of increasing size (in kbp)
windows <- c(500,1000,5000)
bins.for.events <- c(-Inf,0,1,5,10,15,25,50,100,500,Inf)
G2Den.avg.step <- llply(windows, function(window) {
	transform( ddply(G2Map.avg, .(Chr), make.density, window/5, window, .progress = "text"),
						 window = window, events.bin = cut(Events, bins.for.events) )
})
# density in these same windows, but simulating a uniform distribution of events
nsims <- 100
uncertainty.mean <- mean(with(Events, log(End - Start)))
uncertainty.sd <- sd(with(Events, log(End - Start)))
ev.per.chr <- xtabs(~ Chr, Events)
unif.dens.step <- llply(windows, function(window) {
	ldply(1:nsims, function(k) {
		events <- ldply(names(ev.per.chr), function(chr) {
			starts <- round(runif(ev.per.chr[chr], min = min(Events$Start[ Events$Chr == chr ]), max = max(Events$Start[ Events$Chr == chr ])))
			ends <- pmin( starts + round(exp(rnorm(ev.per.chr[chr], mean = uncertainty.mean, sd = uncertainty.sd))) - 1, max(Events$End[ Events$Chr == chr ]) )
			data.frame(Chr = chr, Start = starts, End = ends)
		})
		scalefac.A <- 100 /(Strain.indx %*% Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx)
		scalefac.X <- 100 /(Strain.indx %*% Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx)
		events <- events[ with(events, order(Chr, Start, End)), ]
		map <- ddply(events, .(Chr), integrateCrossovers)
		map$Chr <- as.integer(map$Chr)
		map <- transform(padMap(map), cM = ifelse(Chr < 20, scalefac.A*cumEvents, scalefac.X*cumEvents))
		dens <- ddply(map, .(Chr), make.density, window/5, window)
		dens <- transform(dens, window = window, events.bin = cut(Events, bins.for.events))
		rez <- ddply( ddply(dens, .(window, events.bin), summarize, count = length(events.bin)),
									.(window), transform, prop = count/sum(count) )
		transform(rez, replicate = k)
	}, .parallel = TRUE)
}, .progress = "text")
unif.dens.step <- do.call("rbind", unif.dens.step)
alpha <- 0.000
unif.dens.step.means <- ddply(unif.dens.step, .(window, events.bin), summarize,
															count = mean(count), prop.mean = mean(prop),
															prop.lo = quantile(prop, (alpha/2)), prop.hi = quantile(prop, 1-(alpha/2)))

tmp <- do.call("rbind", G2Den.avg.step)
tmp <- transform(tmp, cM.bin = cut(cM, c(-Inf,0,1e-2,1e-1,1,10)), events.bin = cut(Events, bins.for.events))
# how many windows of size X have less than Y recombination events?
tmp.counts <- ddply( ddply(tmp, .(window, events.bin), summarize, count = length(events.bin)),
										 .(window), transform, prop = count/sum(count) )
# tmp.counts <- rbind(transform(tmp.counts, group = "real data"), transform(unif.dens.step.means, group = "simulations"))
ggplot(tmp.counts) + 
	geom_bar(aes(x = events.bin, y = prop), fill = "grey70", stat = "identity") +
	geom_linerange(data = unif.dens.step.means, aes(x = events.bin, y = prop.mean, ymin = prop.lo, ymax = prop.hi)) +
	geom_point(data = unif.dens.step.means, aes(x = events.bin, y = prop.mean, ymin = prop.lo, ymax = prop.hi), pch = 19, size = 1) +
	geom_line(data = unif.dens.step.means, aes(x = events.bin, y = prop.mean, group = window)) +
	xlab("\nrecombination events per window") + ylab("proportion of windows\n") +
	facet_grid(window ~ ., labeller = function(var,x) paste(x, "kbp")) +
	scale_x_discrete(labels = make.bin.labels(bins.for.events)) +
	theme_bw()

# plot density in windows across genome
qplot(Pos, cM, data=G2Den.avg, geom="line",
      main=paste("G2F1 Avg Density ",as.character(dscale/1000),"Mb",sep="")) +   
      facet_grid(facets=Chr ~.)      
# plot it to file, 1 chromosome per page
limit <- ceiling(max(G2Den.avg$cM))
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_avg.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.avg, Chr==i), geom="line", 
		main=paste("Chromosome", as.character(i), sep=" "))+ylim(c(0,limit)))
    }
dev.off()

## Cox density, for comparison
CoxDen.avg <- ddply(CoxMap.avg, .(Chr), make.density, dstep, dscale)
head(CoxDen.avg)

## plot combined density
CombineDen.avg <- rbind(transform(G2Den.avg, Map="G2F1"), transform(CoxDen.avg, Map="Cox"))
pdf(paste("Density_Cox_G2F1_",as.character(dscale/1000),"Mb.pdf",sep=""))
qplot(Pos, cM, data=CombineDen.avg, geom="line", color=Map,
   main=paste("Avg Density ",as.character(dscale/1000),"Mb",sep="")) +   
   facet_grid(facets=Chr ~., scales="free_y", space="free")
rm(CombineDen.avg)

## find cold regions in sex-averaged G2F1 map using 1d Smith-Waterman-style algorithm (see Materials & Methods)
dscale = 500
mult = 1/100
# compute density
G2Den.avg <- ddply(G2Map.avg, .(Chr), make.density, dscale, dscale)
# find excursions of significantly increased (mult > 1) or decreased (mult < 1) density
G2Den.avg <- ddply(G2Den.avg, .(Chr), make.excursion, mult)
# (do it on bootstrap samples)
G2Den.bs <- llply(G2Map.bs, function(df) {
	df.den <- ddply(df, .(Chr), make.density, dscale, dscale)
	df.den <- ddply(df.den, .(Chr), make.excursion, mult)
	return(df.den)
}, .progress = "text", .parallel = TRUE)

ddply(G2Den.avg[,c(1,6)], .(Chr), numcolwise(max))

## make excursion plot
qplot(Pos, E, data=G2Den.avg, geom="line",
	main=paste("G2F1 Avg Excursion ",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="")) + 
     facet_grid(facets=Chr ~.)

## report the excursion intervals (= cold regions)
(peaks <- ddply(G2Den.avg, .(Chr), get.peaks, 9) )
write.csv(peaks,
	file=paste(paste("G2F1_ColdSpots_avg",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))
# (find cold regions in bootstrap samples)
names(G2Den.bs) <- paste("bootstrap", 1:length(G2Den.bs), sep = ".")
peaks.bs <- ldply(G2Den.bs, function(df) {
	ddply(df, .(Chr), get.peaks, 9)
}, .progress = "text")
	
## now find fold regions in Cox sex-averaged map
CoxDen.avg <- ddply(CoxMap.avg, .(Chr), make.density, dscale, dscale)
CoxDen.avg <- ddply(CoxDen.avg, .(Chr), make.excursion, mult)
ddply(CoxDen.avg[,c(1,6)], .(Chr), numcolwise(max))

## make excursion plot
qplot(Pos, E, data=CoxDen.avg, geom="line") + 
  facet_grid(facets=Chr ~.)  #, scales="free_y", space="free")

## report the intervals
(peaks <- ddply(CoxDen.avg, .(Chr), get.peaks, 34) )
write.csv(peaks,
	file=paste(paste("Cox_ColdSpots_avg",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))

## find hot regions in G2F1 sex-averaged map
dscale = 100
mult = 20
G2Den.avg <- ddply(G2Map.avg, .(Chr), make.density, dscale, dscale) 
G2Den.avg <- ddply(G2Den.avg, .(Chr), make.excursion, mult)

## make excursion plot
qplot(Pos, E, data=G2Den.avg, geom="line",
	main=paste("G2F1 Avg Excursion ",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="")) + 
    facet_grid(facets=Chr ~.)

## report the intervals
(peaks <- ddply(G2Den.avg, .(Chr), get.peaks, 10))
write.csv(peaks,
	file=paste(paste("G2F1_HotSpots_avg",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))

## now find hot regions in Cox sex-averaged map
CoxDen.avg <- ddply(CoxMap.avg, .(Chr), make.density, dscale, dscale)
CoxDen.avg <- ddply(CoxDen.avg, .(Chr), make.excursion, mult)

## make excursion plot
qplot(Pos, E, data=CoxDen.avg, geom="line") + 
  facet_grid(facets=Chr ~.)

## report the intervals
(peaks <- ddply(CoxDen.avg, .(Chr), get.peaks, 40))
write.csv(peaks,
	file=paste(paste("Cox_HotSpots_avg",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))



################################################################
## compute density plots (cM/bp) for sex-specific maps

## total numbers of events by sex, G2F1 vs Cox
total.events.sex <- cbind(
    ddply(subset(G2Map.Sex,Sex=="Female"&Chr<20), .(Chr), numcolwise(max)), 
    ddply(subset(G2Map.Sex,Sex=="Male"), .(Chr), numcolwise(max)),
    ddply(subset(CoxMap.Sex,Sex=="Male"), .(Chr), numcolwise(max)),
    ddply(subset(CoxMap.Sex,Sex=="Female"&Chr<20), .(Chr), numcolwise(max))
    )[,c(1,3,7,12,16)]
names(total.events.sex) <- c("Chr","G2F1.Female.Events","G2F1.Male.Events",	
	"Cox.Female.Events","Cox.Male.Events")
total.events.sex

## prepare to compute density in sliding window: scale = window size; step = offset
# set scale and step size for density computation in kb, 1000=1Mb
dscale = 1000
dstep = dscale/5

## compute density in G2F1 sex-specific maps
G2Den.Sex <- ddply(G2Map.Sex, .(Chr,Sex), make.density, dstep, dscale)
head(G2Den.Sex)
# plot it
qplot(Pos, cM, data=G2Den.Sex, color=Sex, geom="line",
	main=paste("G2F1 Sex Density ",as.character(dscale/1000),"Mb",sep="")) +  
	aes(color=Sex) + scale_colour_manual(value=color.sex) +
    facet_grid(facets=Chr ~ Sex) #, scales="free_y", space="free") 
# plot it to file, 1 chromosome per page
limit <- ceiling(max(G2Den.Sex$cM))
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_Sex.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.Sex, Chr==i), geom="line", 
		main=paste("Chromosome", as.character(i), sep=" "))+ylim(c(0,limit))+
        facet_grid(facets=Sex~.))
    }
dev.off()

## total map length by Sex 		
ddply(ddply(subset(G2Map.Sex,Chr<20),.(Chr, Sex),numcolwise(max))[,c(2,5)],.(Sex),numcolwise(sum))
 
## plot cumulative maps by sex
qplot(Pos, cM, data=G2Map.Sex, geom="line",color=Sex) + 
	aes(color=Sex) + scale_colour_manual(value=color.sex) +
	facet_wrap(facets=~Chr, nrow=5) 
	 
## chr lengths in cM by sex
qplot(male_cM, female_cM, data=subset(ChrLen,Chr<20), colour=Chr) + 
	geom_abline(slope=1,intercept=0,lty=2) +  geom_point(size=3) +
	scale_color_gradientn(colour=rainbow(20), breaks=1:19) + theme_bw()

## plot cumulative sex-speific maps to file, 1 page per chromosome
pdf("CumulativeBySex.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
	print(qplot(Pos, cM, data=subset(G2Map.Sex,Chr==i), 
		geom="line", color=Sex,
		main=paste("Chromosome", as.character(i), sep=" ")) +
		aes(color=Sex) + scale_colour_manual(value=color.sex) )
	}
dev.off()

## now compute density in Cox sex-specific maps
CoxDen.Sex <- ddply(CoxMap.Sex, .(Chr,Sex), make.density, dstep, dscale)
head(CoxDen.Sex)

## plot combined density by sex
CombineDen.Sex <- rbind(transform(G2Den.Sex, Map="G2F1"), transform(CoxDen.Sex, Map="Cox"))
pdf(paste("Density_Cox_G2F1_",as.character(dscale/1000),"Mb_sex_A.pdf",sep=""))
qplot(Pos, cM, data=CombineDen.Sex, geom="line", color=Sex,
   main=paste("Avg Density ",as.character(dscale/1000),"Mb",sep="")) + 
   aes(color=Sex) + scale_colour_manual(value=color.sex) +  
   facet_grid(facets=Chr ~ Map, scales="free_y", space="free")
dev.off()
pdf(paste("Density_Cox_G2F1_",as.character(dscale/1000),"Mb_sex_B.pdf",sep=""))
qplot(Pos, cM, data=CombineDen.Sex, geom="line", color=Map,
   main=paste("Avg Density ",as.character(dscale/1000),"Mb",sep="")) +  
   facet_grid(facets=Chr ~ Sex, scales="free_y", space="free")
dev.off()
# (clean up)
rm(CombineDen.Sex)
 
## find cold regions in G2F1 sex-specific maps 
dscale = 500
mult = 1/100
G2Den.Sex <- ddply(G2Map.Sex, .(Chr,Sex), make.density, dstep, dscale)
G2Den.Sex <- ddply(G2Den.Sex, .(Chr,Sex), make.excursion, mult)
ddply(G2Den.Sex, .(Chr,Sex), numcolwise(max))

## make excursion plot
qplot(Pos, E, data=G2Den.Sex, color=Sex, geom="line") + 
    aes(color=Sex) + scale_colour_manual(value=color.sex) +  
	facet_grid(facets=Chr ~ .)

## report the intervals
(peaks.Female <- ddply(subset(G2Den.Sex,Sex=="Female"), .(Chr), get.peaks, 16) )
write.csv(peaks,
	file=paste(paste("G2F1_ColdSpots_Female",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))
(peaks.Male <- ddply(subset(G2Den.Sex,Sex=="Male"), .(Chr), get.peaks, 16) )
write.csv(peaks,
	file=paste(paste("G2F1_ColdSpots_Male",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))
	
## find cold regions in Cox sex-specific maps
CoxDen.Sex <- ddply(CoxMap.Sex, .(Chr,Sex), make.density, dstep, dscale)
CoxDen.Sex <- ddply(CoxDen.Sex, .(Chr,Sex), make.excursion, mult)
ddply(CoxDen.Sex[,c(1,6)], .(Chr,Sex), numcolwise(max))

## make excursion plot
qplot(Pos, E, data=CoxDen.Sex, color=Sex, geom="line") +     
	aes(color=Sex) + scale_colour_manual(value=color.sex) +  
    facet_grid(facets=Chr ~.)  #, scales="free_y", space="free")

## report the intervals
(peaks.Female <- ddply(subset(CoxDen.Sex, Sex=="Female"), .(Chr), get.peaks, 60) )
write.csv(peaks,
	file=paste(paste("Cox_ColdSpots_Female",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))
(peaks.Male <- ddply(subset(CoxDen.Sex, Sex=="Male"), .(Chr), get.peaks, 60) )
write.csv(peaks,
	file=paste(paste("Cox_ColdSpots_Male",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))

## find hot regions in G2F1 sex-specific maps
dscale = 100
mult = 20
G2Den.Sex <- ddply(G2Map.Sex, .(Chr,Sex), make.density, dscale, dscale)
G2Den.Sex <- ddply(G2Den.Sex, .(Chr,Sex), make.excursion, mult)

## make excursion plot
qplot(Pos, E, data=G2Den.Sex, color=Sex, geom="line") + 
  aes(color=Sex) + scale_colour_manual(value=color.sex) +  
  facet_grid(facets=Chr ~.)

## report the intervals
(peaks.Female <- ddply(subset(G2Den.Sex,Sex=="Female"), .(Chr), get.peaks, 8))
write.csv(peaks,
	file=paste(paste("G2F1_HotSpots_Female",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))
(peaks.Male <- ddply(subset(G2Den.Sex,Sex=="Male"), .(Chr), get.peaks, 11))
write.csv(peaks,
	file=paste(paste("G2F1_HotSpots_Male",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))

## now find hot regions in Cox sex-specific maps
CoxDen.Sex <- ddply(CoxMap.Sex, .(Chr,Sex), make.density, dscale, dscale)
CoxDen.Sex <- ddply(CoxDen.Sex, .(Chr,Sex), make.excursion, mult)

## make excursion plot
qplot(Pos, E, data=CoxDen.Sex, color=Sex, geom="line") + 
  aes(color=Sex) + scale_colour_manual(value=color.sex) +  
  facet_grid(facets=Chr ~.) 

## report the intervals
(peaks.Female <- ddply(subset(CoxDen.Sex, Sex=="Female"), .(Chr), get.peaks, 28))
write.csv(peaks,
	file=paste(paste("Cox_HotSpots_Female",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))
(peaks.Male <- ddply(subset(CoxDen.Sex, Sex=="Male"), .(Chr), get.peaks, 38))
write.csv(peaks,
	file=paste(paste("Cox_HotSpots_Male",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))



####################################################
## Compute strain-specific cumulative maps in G2F1 population
# a recombination between strains A,B counts in the tally for both A and B

G2Map.Strain <- rbind(
  transform(padMap(ddply(subset(Events,S129),.(Chr),integrateCrossovers)),Strain="S129"),  
   transform(padMap(ddply(subset(Events,AJ),.(Chr),integrateCrossovers)),Strain="AJ"),  
   transform(padMap(ddply(subset(Events,B6),.(Chr),integrateCrossovers)),Strain="B6"),  
   transform(padMap(ddply(subset(Events,CAST),.(Chr),integrateCrossovers)),Strain="CAST"),  
   transform(padMap(ddply(subset(Events,NOD),.(Chr),integrateCrossovers)),Strain="NOD"),  
   transform(padMap(ddply(subset(Events,NZO),.(Chr),integrateCrossovers)),Strain="NZO"),  
   transform(padMap(ddply(subset(Events,PWK),.(Chr),integrateCrossovers)),Strain="PWK"),  
   transform(padMap(ddply(subset(Events,WSB),.(Chr),integrateCrossovers)),Strain="WSB"))  
# (peek to check result)
head(G2Map.Strain)

## scale the cumulative maps to cM units
scalefac.A <- 100 /Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx
scalefac.X <- 100 /Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx
G2Map.Strain <- transform(G2Map.Strain, 
	cM = ifelse(Chr<20, cumEvents*scalefac.A[Strain], cumEvents*scalefac.X[Strain] ))
# (peek to check result)
head(G2Map.Strain)

## check that the average chr length across strains is the same as avg chr length
#  X will be a little off due to unequal weights
ddply(ddply(G2Map.Strain,.(Chr, Strain),numcolwise(max))[,c(1,5)],.(Chr),numcolwise(mean))[,2] /
	(ddply(G2Map.avg,.(Chr),numcolwise(max))[,4])

## check total length of map by strain	
ddply(ddply(G2Map.Strain,.(Chr, Strain),numcolwise(max))[,c(2,5)],.(Strain),numcolwise(sum))
#  Strain       cM
#1   S129 1383.244
#2     AJ 1359.956
#3     B6 1385.954
#4   CAST 1364.612
#5    NOD 1386.709
#6    NZO 1376.966
#7    PWK 1336.163
#8    WSB 1358.674

## update chromosome lengths
ChrLen <- transform(ChrLen,
	Prdm9_A_cM = ddply(subset(G2Map.Prdm9,Allele=="Prdm9_A"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_B_cM = ddply(subset(G2Map.Prdm9,Allele=="Prdm9_B"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_C_cM = ddply(subset(G2Map.Prdm9,Allele=="Prdm9_C"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_D_cM = ddply(subset(G2Map.Prdm9,Allele=="Prdm9_D"),.(Chr),	
		numcolwise(max))[,"cM"]
)
ChrLen 

## plot physical x genetic size by strain
qplot(Pos/1e6,cM,
	data=subset(ddply(G2Map.Strain,.(Chr, Strain),numcolwise(max)),Chr<20) ) +
	aes(colour=Strain)+scale_colour_manual(values=CC.colors.alt)	+ 
	geom_point(aes(Pos/1e6,cM),data=subset(ddply(G2Map.Strain,.(Chr, Strain),
		numcolwise(max)),Chr==20)) + xlab("Size in Mb") + ylab("Size in cM") 

## plot cumulative maps by strain, 1 page per chromosome, to a file
pdf("CumulativeByStrain.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
	print(qplot(Pos, cM, data=subset(G2Map.Strain,Chr==i), geom="line", 
			color=Strain,main=paste("Chromosome", as.character(i), sep=" ")) +
			aes(colour=Strain)+scale_colour_manual(values=CC.colors.alt)  )
}	
dev.off()
pdf("CumulativeByStrain2.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
	print(qplot(Pos, cM, data=subset(G2Map.Strain,Chr==i), geom="line", 
			color=Strain,main=paste("Chromosome", as.character(i), sep=" ")) +
			aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt) + 
			facet_wrap(facets=~Strain, nrow=4)
 )
}
dev.off()

## compute recombination density by strain
dscale = 100
dstep = dscale/5
G2Den.Strain <- ddply(G2Map.Strain, .(Chr, Strain), make.density, dstep, dscale)
head(G2Den.Strain)

## plot density to file, 1 page per chromosome
limit = ceiling(max(G2Den.Strain$cM))    	
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_strain.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.Strain, Chr==i), geom="line", 
		main=paste("Chromosome", as.character(i), sep=" ")) + 
    	aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt)+
    	ylim(c(0,limit)) + facet_grid(facets=Strain ~.) )
    }
dev.off()

## Figure S8: a PWK-specific hotspot on distal chr15
ggplot(subset(G2Den.Strain, Chr == 15 & Pos > 95e6 & Pos < 105e6)) +
	geom_segment(aes(x = Pos/1e6, xend = Pos/1e6, y = 0, yend = cM, colour = Strain)) + facet_grid(Strain~Chr) +
	scale_colour_manual(values = CC.colors.alt) +
	xlab("\ncoordinate (Mbp)") + ylab("\ncM per 20 kbp\n")

## find cold regions in strain-specific G2F1 maps
dscale = 500
mult = 1/100
G2Den.Strain <- ddply(G2Map.Strain, .(Chr,Strain), make.density, dscale, dscale)
G2Den.Strain <- ddply(G2Den.Strain, .(Chr,Strain), make.excursion, mult)
ddply(G2Den.Strain, .(Chr,Strain), numcolwise(max))

## make excursion plot
qplot(Pos, E, data=G2Den.Strain, color=Strain, geom="line") + 
    aes(color=Strain) + scale_colour_manual(value=CC.colors.alt) +  
	facet_grid(facets=Chr ~ .)  

## report the intervals per strain
peaks <- NULL
for(name in StrainNames){
	peaks <- rbind(peaks,
		transform(ddply(subset(G2Den.Strain,Strain==name), .(Chr), get.peaks, 10),
		Strain=name) )
	}
write.csv(peaks,
	file=paste(paste("G2F1_ColdSpots_Strain",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))

## find hot regions in strain-specific G2F1 maps
dscale = 100
mult = 20
G2Den.Strain <- ddply(G2Map.Strain, .(Chr,Strain), make.density, dscale, dscale, .progress = "text")
G2Den.Strain <- ddply(G2Den.Strain, .(Chr,Strain), make.excursion, mult, .progress = "text")
ddply(G2Den.Strain, .(Chr,Strain), numcolwise(max))
# (now shuffle strain labels and compute background distribution??)
## now do strain-specific maps on bootstrap samples: in each sample, permute strain labels ONCE
# use plyr's llply so I get a progress bar
library(doMC)
registerDoMC(cores = 8)
n.bs.reps <- 100
bs.quantile <- 0.99
peaks.bs <- ldply(1:n.bs.reps, function(i) {
	
	# local copy of Events table
	df <- Events
	
	# first shuffle strains at junction
	df$Dist_Founder <- sample(df$Dist_Founder)
	df$Prox_Founder <- sample(df$Prox_Founder)
	# then re-create indicator variables
	df <- transform(df,
									S129 = (Dist_Founder=="S129" | Prox_Founder=="S129"),
									AJ = (Dist_Founder=="AJ" | Prox_Founder=="AJ"),
									B6 = (Dist_Founder=="B6" | Prox_Founder=="B6"),
									CAST = (Dist_Founder=="CAST" | Prox_Founder=="CAST"),
									NOD = (Dist_Founder=="NOD" | Prox_Founder=="NOD"), 
									NZO = (Dist_Founder=="NZO" | Prox_Founder=="NZO"),
									PWK = (Dist_Founder=="PWK" | Prox_Founder=="PWK"),
									WSB = (Dist_Founder=="WSB" | Prox_Founder=="WSB") )
	
	# make strain-specific cumulative maps
	map.strain <- rbind(
		transform(padMap(ddply(subset(df,S129),.(Chr),integrateCrossovers)),Strain="S129"),  
		transform(padMap(ddply(subset(df,AJ),.(Chr),integrateCrossovers)),Strain="AJ"),  
		transform(padMap(ddply(subset(df,B6),.(Chr),integrateCrossovers)),Strain="B6"),  
		transform(padMap(ddply(subset(df,CAST),.(Chr),integrateCrossovers)),Strain="CAST"),  
		transform(padMap(ddply(subset(df,NOD),.(Chr),integrateCrossovers)),Strain="NOD"),  
		transform(padMap(ddply(subset(df,NZO),.(Chr),integrateCrossovers)),Strain="NZO"),  
		transform(padMap(ddply(subset(df,PWK),.(Chr),integrateCrossovers)),Strain="PWK"),  
		transform(padMap(ddply(subset(df,WSB),.(Chr),integrateCrossovers)),Strain="WSB"))
	
	# get scaling factors...
	sf <- get.scaling.factors(df, StrainNames, MeiosisNames)
	scalefac.A <- 100 /sf$Strain.Pos %*% sf$Pos.Meiosis.A %*% sf$Meiosis.indx
	scalefac.X <- 100 /sf$Strain.Pos %*% sf$Pos.Meiosis.X %*% sf$Meiosis.indx
	# ...and rescale
	map.strain <- transform(map.strain, 
													cM = ifelse(Chr<20, cumEvents*scalefac.A[Strain], cumEvents*scalefac.X[Strain] ))
	
	# do density calculation
	den.strain <- ddply(map.strain, .(Chr,Strain), make.density, dscale, dscale)
	den.strain <- ddply(den.strain, .(Chr,Strain), make.excursion, mult)
	
	# ddply(den.strain, .(Chr,Strain), summarize, E.max = quantile(E, bs.quantile))
	ddply(den.strain, .(Chr,Strain), get.peaks, 12)
	
}, .progress = "text", .parallel = TRUE)
# G2Den.Strain.bg <- llply(G2Map.Strain.bs, function(df) {
# 	df.den <- ddply(df, .(Chr,Strain), make.density, dscale, dscale)
# 	df.den <- ddply(df.den, .(Chr,Strain), make.excursion, mult)
# }, .progress = "text", .parallel = TRUE)

## make excursion plot
qplot(Pos, E, data=G2Den.Strain, color=Strain, geom="line") + 
    aes(color=Strain) + scale_colour_manual(value=CC.colors.alt) +  
	facet_grid(facets=Chr ~ .)  

## report the intervals per strain
peaks <- NULL
for(name in StrainNames){
	strain.peaks <- ddply(subset(G2Den.Strain,Strain==name), .(Chr), get.peaks, 12)
	if (nrow(strain.peaks)) {
		peaks <- rbind(peaks, transform(strain.peaks, Strain=name))
	}
}
write.csv(peaks,
	file=paste(paste("G2F1_HotSpots_Strain",dscale,"kbp",
	"x",as.character(mult),sep="_"),".csv",sep=""))
## plot real peaks versus permutation E-scores
peaks.all <- rbind( transform(peaks, group = "real data"), transform(peaks.bs, group = "permutations") )
peaks.all$Strain <- factor(peaks.all$Strain, levels = StrainNames[ c(2,3,1,5,6,4,7,8) ])
ggplot(peaks.all) +
	geom_boxplot(aes(y = E, x = Strain, fill = Strain, alpha = (group != "permutations"))) +
	scale_fill_manual(values = CC.colors.alt) + scale_colour_manual(values = CC.colors.alt) +
	scale_alpha_manual(values = c(0.3,1.0)) + 
	ylab("segment E-score\n") + xlab("") + ylim(c(0,50)) +
	guides(fill = FALSE, alpha = FALSE) +
	theme_bw()

## (save distribution in bootstrap samples)
# names(G2Den.Strain.bs) <- paste("bootstrap", 1:length(G2Den.Strain.bs), sep = ".")
# peaks.strain.bs <- dlply(G2Den.Strain.bs, function(df) {
# 	peaks <- NULL
# 	for(name in StrainNames){
# 		peaks <- rbind(peaks,
# 									 transform(ddply(subset(G2Den.Strain,Strain==name), .(Chr), get.peaks, 12),
# 									 					Strain=name) )
# 	}
# 	return(peaks)
# }, .progress = "text", .parallel = TRUE)


##########################################################
## Compute cumulative maps in G2F1 population, excluding breakpoints involving one strain at a time
# in this analysis, a recombination event counts in the tally for not-A if strain A is not on either side of breakpoint
# NB: this is mostly a sanity check on our inference of founder haplotypes

G2Map.notStrain <- rbind(
  transform(padMap(ddply(subset(Events,!S129),.(Chr),integrateCrossovers)),Strain="S129"),  
  transform(padMap(ddply(subset(Events,!AJ),.(Chr),integrateCrossovers)),Strain="AJ"),  
  transform(padMap(ddply(subset(Events,!B6),.(Chr),integrateCrossovers)),Strain="B6"),  
  transform(padMap(ddply(subset(Events,!CAST),.(Chr),integrateCrossovers)),Strain="CAST"),  
  transform(padMap(ddply(subset(Events,!NOD),.(Chr),integrateCrossovers)),Strain="NOD"),  
  transform(padMap(ddply(subset(Events,!NZO),.(Chr),integrateCrossovers)),Strain="NZO"),  
  transform(padMap(ddply(subset(Events,!PWK),.(Chr),integrateCrossovers)),Strain="PWK"),  
  transform(padMap(ddply(subset(Events,!WSB),.(Chr),integrateCrossovers)),Strain="WSB"))  

## scale cumulative maps to cM units
scalefac.A <- 100 /(7*n.funnels - Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx)
scalefac.X <- 100 /(3.25*n.funnels - Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx)
G2Map.notStrain <- transform(G2Map.notStrain, 
	cM = ifelse(Chr<20, cumEvents*scalefac.A[Strain], cumEvents*scalefac.X[Strain] ))
# (peek to check result)
head(G2Map.notStrain)
 
## double-check the scaling
ddply(ddply(G2Map.notStrain,.(Chr, Strain), numcolwise(max)),.(Chr),numcolwise(mean))$cM /
  		ddply(G2Map.avg, .(Chr), numcolwise(max))$cM
  		
## total map length by excluded strain 		
ddply(ddply(G2Map.notStrain,.(Chr, Strain),numcolwise(max))[,c(2,5)],.(Strain),numcolwise(sum))
#  Strain       cM
#1   S129 1364.420
#2     AJ 1372.370
#3     B6 1363.585
#4   CAST 1370.557
#5    NOD 1363.268
#6    NZO 1366.267
#7    PWK 1379.986
#8    WSB 1372.615

## plot the cumulative maps per excluded strain, 1 page per chromosome, to file
pdf("CumulativeByNotStrain.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
	print(qplot(Pos, cM, data=subset(rbind(
		transform(G2Map.Strain, Include="Strain"),
		transform(G2Map.notStrain, Include="notStrain")), Chr==i), 
		geom="line", lty=Include,color=Strain,
		main=paste("Chromosome", as.character(i), sep=" ")) +
		aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt) + 
		facet_wrap(facets=~Strain, nrow=4)  )
}
dev.off()



####################################################
## Compute sex- AND strain-specific cumulative maps in G2F1 population

Events.male <- Events[apply(outer(Events[,"Type"],c("MGP","PGP","P"),"=="),1,any),]
Events.female <- Events[apply(outer(Events[,"Type"],c("MGM","PGM","M"),"=="),1,any),]

## prepare scaling factors
scalefac.A <- 100 /Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx.sex
scalefac.X <- 100 /Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx.sex[,2]

## compute maps per sex-strain combination
G2Map.SexStrain <- NULL
for(i in 1:8){
	G2Map.SexStrain <- rbind(G2Map.SexStrain,
	transform(padMap(ddply(Events.male[Events.male[,32+i],],.(Chr),integrateCrossovers)),	
   		Sex="Male",Strain=StrainNames[i],cM=cumEvents*scalefac.A[i,1]))
	}
for(i in 1:8){
	G2Map.SexStrain <- rbind(G2Map.SexStrain,
	transform(padMap(ddply(Events.female[Events.female[,32+i],],.(Chr),integrateCrossovers)),	
   		Sex="Female",Strain=StrainNames[i],
   		cM=ifelse(Chr<20,cumEvents*scalefac.A[i,2],cumEvents*scalefac.X[i])))
	}
# (peek to check result)
head(G2Map.SexStrain)

## total map lengths by sex and strain
ddply(ddply(subset(G2Map.SexStrain,Chr<20),.(Chr, Sex, Strain),	
	numcolwise(max))[,c(2,3,6)],.(Sex,Strain),numcolwise(sum))

## double-check scaling
ddply(ddply(G2Map.SexStrain,.(Chr, Sex, Strain),numcolwise(max)),.(Chr),numcolwise(mean))$cM /
  ddply(G2Map.avg,.(Chr),numcolwise(max))$cM
  
## update chromosome lengths
ChrLen <- transform(ChrLen,
	S129_Male_cM = c(ddply(subset(G2Map.SexStrain,Strain=="S129"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	AJ_Male_cM = c(ddply(subset(G2Map.SexStrain,Strain=="AJ"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	B6_Male_cM = c(ddply(subset(G2Map.SexStrain,Strain=="B6"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	CAST_Male_cM = c(ddply(subset(G2Map.SexStrain,Strain=="CAST"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	NOD_Male_cM = c(ddply(subset(G2Map.SexStrain,Strain=="NOD"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	NZO_Male_cM = c(ddply(subset(G2Map.SexStrain,Strain=="NZO"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	PWK_Male_cM = c(ddply(subset(G2Map.SexStrain,Strain=="PWK"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),
	WSB_Male_cM = c(ddply(subset(G2Map.SexStrain,Strain=="WSB"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),
	S129_Female_cM = ddply(subset(G2Map.SexStrain,Strain=="S129"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	AJ_Female_cM = ddply(subset(G2Map.SexStrain,Strain=="AJ"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	B6_Female_cM = ddply(subset(G2Map.SexStrain,Strain=="B6"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	CAST_Female_cM = ddply(subset(G2Map.SexStrain,Strain=="CAST"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	NOD_Female_cM = ddply(subset(G2Map.SexStrain,Strain=="NOD"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	NZO_Female_cM = ddply(subset(G2Map.SexStrain,Strain=="NZO"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	PWK_Female_cM = ddply(subset(G2Map.SexStrain,Strain=="PWK"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],
	WSB_Female_cM = ddply(subset(G2Map.SexStrain,Strain=="WSB"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"]	
)
ChrLen

## plot cumulative maps by sex and strain, 1 page per chromosome, to file
pdf("CumulativeBySexStrain.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
	print(qplot(Pos, cM, data=subset(G2Map.SexStrain,Chr==i), geom="line", 
			color=Strain,main=paste("Chromosome", as.character(i), sep=" ")) +
			aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt) +
			facet_grid(Sex~.))
	}
dev.off()		
pdf("CumulativeByStrainSex.pdf", width=9.5, height=7, onefile= TRUE)			
for(i in 1:20){
	print(qplot(Pos, cM, data=subset(G2Map.SexStrain,Chr==i), geom="line", linetype=Sex,
			color=Strain,main=paste("Chromosome", as.character(i), sep=" ")) +
			aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt) +
			facet_wrap(facets=~Strain, nrow=4))
}		
dev.off()

## compute recombination density by sex and strain
dscale = 100
dstep = dscale/5
G2Den.SexStrain <- ddply(G2Map.SexStrain, .(Chr, Sex, Strain), make.density, dstep, dscale)
head(G2Den.SexStrain)

## plot density to file, 1 page per chromosome
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_SexStrain.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexStrain, Chr==i), geom="line",
		main=paste("Chromosome", as.character(i), sep=" ")) + 
    	aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt)+
    	facet_grid(Strain~Sex))
    }
dev.off()
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_StrainSex.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexStrain, Chr==i), geom="line", linetype=Sex,
		main=paste("Chromosome", as.character(i), sep=" ")) + 
    	aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt)+
    	facet_wrap(facets=~Strain, nrow=4) )
    }
dev.off()
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_SexStrain2.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexStrain, Chr==i), geom="line",
		main=paste("Chromosome", as.character(i), sep=" ")) + 
    	aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt)+
    	facet_grid(Sex~.))
}
dev.off()

## find sex- and strain-specific cold regions in G2F1 map
dscale = 500
mult = 1/100
G2Den.SexStrain <- ddply(G2Map.SexStrain, .(Chr,Sex,Strain), make.density, dscale, dscale)
G2Den.SexStrain <- ddply(G2Den.SexStrain, .(Chr,Sex,Strain), make.excursion, mult)
ddply(G2Den.SexStrain, .(Chr,Sex, Strain), numcolwise(max))

## make excursion plot
qplot(Pos, E, data=G2Den.SexStrain, color=Strain, geom="line") + 
    aes(color=Strain) + scale_colour_manual(value=CC.colors.alt) +  
	facet_grid(facets=Chr ~ Sex)  

## report the intervals per strain for MALES
peaks.Male <- NULL
for(name in StrainNames){
	peaks.Male <- rbind(peaks,
		transform(ddply(subset(G2Den.SexStrain,Sex=="Male"&Strain==name), 
		.(Chr), get.peaks, 80),
		Strain=name) )
	}
write.csv(peaks.Male,
	file=paste(paste("G2F1_Male_ColdSpots_SexStrain",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))

## report the intervals per strain for FEMALES
peaks.Female <- NULL
for(name in StrainNames){
	peaks.Female <- rbind(peaks,
		transform(ddply(subset(G2Den.SexStrain,Sex=="Female"&Strain==name), 
		.(Chr), get.peaks, 10),
		Strain=name) )
}
write.csv(peaks,
	file=paste(paste("G2F1_Female_ColdSpots_SexStrain",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))

## find sex- and strain-specific hot regions in G2F1 map
dscale = 100
mult = 20
G2Den.SexStrain <- ddply(G2Map.SexStrain, .(Chr,Strain,Sex), make.density, dscale, dscale)
G2Den.SexStrain <- ddply(G2Den.SexStrain, .(Chr,Strain,Sex), make.excursion, mult)
ddply(G2Den.SexStrain, .(Chr,Strain,Sex), numcolwise(max))

## make excursion plot
qplot(Pos, E, data=G2Den.SexStrain, color=Strain, geom="line") + 
    aes(color=Strain) + scale_colour_manual(value=CC.colors.alt) +  
	facet_grid(facets=Chr ~ Sex)  

## report the intervals per strain for MALES
peaks.Male <- NULL
for(name in StrainNames){
	peaks <- rbind(peaks,
		transform(ddply(subset(G2Den.SexStrain,Sex=="Male"&Strain==name), 
		.(Chr), get.peaks, 10),
		Strain=name, Sex="Male") )
}	
write.csv(peaks,
	file=paste(paste("G2F1_Male_HotSpots_SexStrain",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))

## report the intervals per strain for FEMALES
peaks.Female <- NULL
for(name in StrainNames){
	peaks <- rbind(peaks,
		transform(ddply(subset(G2Den.SexStrain,Sex=="Female"&Strain==name), 
		.(Chr), get.peaks, 10),
		Strain=name, Sex="Female") )
}	
write.csv(peaks,
	file=paste(paste("G2F1_Female_HotSpots_SexStrain",as.character(dscale/1000),"Mb",
	"x",as.character(mult),sep="_"),".csv",sep=""))



###########################################
## Compute sex- AND strain- AND generation-specific maps in G2F1 population 

## Figures 2,S6: assign X-chromosome status to G1-male events, given funnel order (see Materials & Methods; Figure 1; and Supplementary Text)
funnel.orders <- matrix(unlist(strsplit(levels(Events$Funnel_Order), "")),
												ncol = 8, byrow = TRUE)
rownames(funnel.orders) <- levels(Events$Funnel_Order)
X.idx <- setNames(c(3,7), c("MGP","PGP"))
Events$X.strain <- apply(Events, 1, function(row) {
	if (as.character(row["Type"]) %in% c("MGP","PGP")) {
		f.pos <- X.idx[ as.character(row["Type"]) ]
		funnel.orders[ row["Funnel_Order"], f.pos ]
	} else {
		NA
	}
})
Events <- transform( Events, X.strain = ifelse(X.strain == "F", "CAST",
																				ifelse(X.strain == "G", "MUS",
																 					ifelse(!is.na(X.strain), "DOM", NA))) )

## prepare scaling factors for converting cumulative maps to cM units
Meiosis.indx.sexgen <- matrix(
	c(0,1,0,1,0,0,
	  1,0,1,0,0,0,
	  0,0,0,0,0,1,
	  0,0,0,0,1,0),ncol=4)
rownames(Meiosis.indx.sexgen) <- MeiosisNames
colnames(Meiosis.indx.sexgen) <- c("G1_Male","G1_Female","G2_Male","G2_Female")
scalefac.A <- 100 /Strain.Pos %*% Pos.Meiosis.A %*% Meiosis.indx.sexgen
scalefac.X <- 100 /Strain.Pos %*% Pos.Meiosis.X %*% Meiosis.indx.sexgen[,c(2,4)]

## Figure 2: scaling factors for map length in G1, by X-chromosome origin
G1X.Pos <- c(1,0,1,0,1,0,1,0)
scalefac.G1X <- Strain.Pos %*% G1X.Pos
scalefac.G1X.subspec <- c( sum(scalefac.G1X[4,]), sum(scalefac.G1X[ c(1:3,5,6,8), ]), sum(scalefac.G1X[7,]) )
names(scalefac.G1X.subspec) <- c("CAST","DOM","MUS")
Meiosis.indx.G1X <- matrix(c(1,1,0,0,0,0,0,0,
														 0,0,1,1,0,0,0,0,
														 0,0,0,0,1,1,0,0,
														 0,0,0,0,0,0,1,1), nrow = 4, byrow = TRUE)
scalefac.G1A <- Strain.Pos %*% t(Meiosis.indx.G1X) %*% t(Meiosis.indx.sexgen)
scalefac.G1A.subspec <- rbind( CAST = colSums(t(as.matrix(scalefac.G1A[4,]))),
															 DOM = colSums(scalefac.G1A[ c(1:3,5,6,8), ]),
															 MUS = colSums(t(as.matrix(scalefac.G1A[7,]))) )

## double-check scaling
ddply(ddply(G2Map.SexStrainGen,.(Chr, Sex, Strain, Gen),numcolwise(max)),
			.(Chr),numcolwise(mean))$cM / ddply(G2Map.avg,.(Chr),numcolwise(max))$cM
ddply(ddply(subset(G2Map.SexStrainGen,Gen=="G1"),.(Chr,Sex,Strain),numcolwise(max)),
			.(Chr),numcolwise(mean))$cM / ddply(G2Map.avg,.(Chr),numcolwise(max))$cM 
ddply(ddply(subset(G2Map.SexStrainGen,Gen=="G2"),.(Chr,Sex,Strain),numcolwise(max)),
			.(Chr),numcolwise(mean))$cM / ddply(G2Map.avg,.(Chr),numcolwise(max))$cM  

## update chromosome lengths
ChrLen <- transform(ChrLen,
										S129_Male_Gen1_cM = c(ddply(subset(G2Map.SexStrainGen,
																											 Strain=="S129"&Sex=="Male"&Gen=="G1"),.(Chr),	
																								numcolwise(max))[,"cM"],NA),	
										AJ_Male_Gen1_cM = c(ddply(subset(G2Map.SexStrainGen,
																										 Strain=="AJ"&Sex=="Male"&Gen=="G1"),.(Chr),	
																							numcolwise(max))[,"cM"],NA),	
										B6_Male_Gen1_cM = c(ddply(subset(G2Map.SexStrainGen,
																										 Strain=="B6"&Sex=="Male"&Gen=="G1"),.(Chr),	
																							numcolwise(max))[,"cM"],NA),	
										CAST_Male_Gen1_cM = c(ddply(subset(G2Map.SexStrainGen,
																											 Strain=="CAST"&Sex=="Male"&Gen=="G1"),.(Chr),	
																								numcolwise(max))[,"cM"],NA),	
										NOD_Male_Gen1_cM = c(ddply(subset(G2Map.SexStrainGen,
																											Strain=="NOD"&Sex=="Male"&Gen=="G1"),.(Chr),	
																							 numcolwise(max))[,"cM"],NA),	
										NZO_Male_Gen1_cM = c(ddply(subset(G2Map.SexStrainGen,
																											Strain=="NZO"&Sex=="Male"&Gen=="G1"),.(Chr),	
																							 numcolwise(max))[,"cM"],NA),	
										PWK_Male_Gen1_cM = c(ddply(subset(G2Map.SexStrainGen,
																											Strain=="PWK"&Sex=="Male"&Gen=="G1"),.(Chr),	
																							 numcolwise(max))[,"cM"],NA),
										WSB_Male_Gen1_cM = c(ddply(subset(G2Map.SexStrainGen,
																											Strain=="WSB"&Sex=="Male"&Gen=="G1"),.(Chr),	
																							 numcolwise(max))[,"cM"],NA),
										S129_Female_Gen1_cM = ddply(subset(G2Map.SexStrainGen,
																											 Strain=="S129"&Sex=="Female"&Gen=="G1"),.(Chr),	
																								numcolwise(max))[,"cM"],	
										AJ_Female_Gen1_cM = ddply(subset(G2Map.SexStrainGen,
																										 Strain=="AJ"&Sex=="Female"&Gen=="G1"),.(Chr),	
																							numcolwise(max))[,"cM"],	
										B6_Female_Gen1_cM = ddply(subset(G2Map.SexStrainGen,
																										 Strain=="B6"&Sex=="Female"&Gen=="G1"),.(Chr),	
																							numcolwise(max))[,"cM"],	
										CAST_Female_Gen1_cM = ddply(subset(G2Map.SexStrainGen,
																											 Strain=="CAST"&Sex=="Female"&Gen=="G1"),.(Chr),	
																								numcolwise(max))[,"cM"],	
										NOD_Female_Gen1_cM = ddply(subset(G2Map.SexStrainGen,
																											Strain=="NOD"&Sex=="Female"&Gen=="G1"),.(Chr),	
																							 numcolwise(max))[,"cM"],	
										NZO_Female_Gen1_cM = ddply(subset(G2Map.SexStrainGen,
																											Strain=="NZO"&Sex=="Female"&Gen=="G1"),.(Chr),	
																							 numcolwise(max))[,"cM"],	
										PWK_Female_Gen1_cM = ddply(subset(G2Map.SexStrainGen,
																											Strain=="PWK"&Sex=="Female"&Gen=="G1"),.(Chr),	
																							 numcolwise(max))[,"cM"],
										WSB_Female_Gen1_cM = ddply(subset(G2Map.SexStrainGen,
																											Strain=="WSB"&Sex=="Female"&Gen=="G1"),.(Chr),	
																							 numcolwise(max))[,"cM"],
										S129_Male_Gen2_cM = c(ddply(subset(G2Map.SexStrainGen,
																											 Strain=="S129"&Sex=="Male"&Gen=="G2"),.(Chr),	
																								numcolwise(max))[,"cM"],NA),	
										AJ_Male_Gen2_cM = c(ddply(subset(G2Map.SexStrainGen,
																										 Strain=="AJ"&Sex=="Male"&Gen=="G2"),.(Chr),	
																							numcolwise(max))[,"cM"],NA),	
										B6_Male_Gen2_cM = c(ddply(subset(G2Map.SexStrainGen,
																										 Strain=="B6"&Sex=="Male"&Gen=="G2"),.(Chr),	
																							numcolwise(max))[,"cM"],NA),	
										CAST_Male_Gen2_cM = c(ddply(subset(G2Map.SexStrainGen,
																											 Strain=="CAST"&Sex=="Male"&Gen=="G2"),.(Chr),	
																								numcolwise(max))[,"cM"],NA),	
										NOD_Male_Gen2_cM = c(ddply(subset(G2Map.SexStrainGen,
																											Strain=="NOD"&Sex=="Male"&Gen=="G2"),.(Chr),	
																							 numcolwise(max))[,"cM"],NA),	
										NZO_Male_Gen2_cM = c(ddply(subset(G2Map.SexStrainGen,
																											Strain=="NZO"&Sex=="Male"&Gen=="G2"),.(Chr),	
																							 numcolwise(max))[,"cM"],NA),	
										PWK_Male_Gen2_cM = c(ddply(subset(G2Map.SexStrainGen,
																											Strain=="PWK"&Sex=="Male"&Gen=="G2"),.(Chr),	
																							 numcolwise(max))[,"cM"],NA),
										WSB_Male_Gen2_cM = c(ddply(subset(G2Map.SexStrainGen,
																											Strain=="WSB"&Sex=="Male"&Gen=="G2"),.(Chr),	
																							 numcolwise(max))[,"cM"],NA),
										S129_Female_Gen2_cM = ddply(subset(G2Map.SexStrainGen,
																											 Strain=="S129"&Sex=="Female"&Gen=="G2"),.(Chr),	
																								numcolwise(max))[,"cM"],	
										AJ_Female_Gen2_cM = ddply(subset(G2Map.SexStrainGen,
																										 Strain=="AJ"&Sex=="Female"&Gen=="G2"),.(Chr),	
																							numcolwise(max))[,"cM"],	
										B6_Female_Gen2_cM = ddply(subset(G2Map.SexStrainGen,
																										 Strain=="B6"&Sex=="Female"&Gen=="G2"),.(Chr),	
																							numcolwise(max))[,"cM"],	
										CAST_Female_Gen2_cM = ddply(subset(G2Map.SexStrainGen,
																											 Strain=="CAST"&Sex=="Female"&Gen=="G2"),.(Chr),	
																								numcolwise(max))[,"cM"],	
										NOD_Female_Gen2_cM = ddply(subset(G2Map.SexStrainGen,
																											Strain=="NOD"&Sex=="Female"&Gen=="G2"),.(Chr),	
																							 numcolwise(max))[,"cM"],	
										NZO_Female_Gen2_cM = ddply(subset(G2Map.SexStrainGen,
																											Strain=="NZO"&Sex=="Female"&Gen=="G2"),.(Chr),	
																							 numcolwise(max))[,"cM"],	
										PWK_Female_Gen2_cM = ddply(subset(G2Map.SexStrainGen,
																											Strain=="PWK"&Sex=="Female"&Gen=="G2"),.(Chr),	
																							 numcolwise(max))[,"cM"],
										WSB_Female_Gen2_cM = ddply(subset(G2Map.SexStrainGen,
																											Strain=="WSB"&Sex=="Female"&Gen=="G2"),.(Chr),	
																							 numcolwise(max))[,"cM"]		
)
ChrLen

#### reproduce some tables and figures from main text ####
## Table 3: construct bootstrapped maps by sex and generation
# how many replicates?
nreps <- 100
G2Map.bootstrap <- ldply(1:nreps,
												 function(j) {
												 	
												 	## draw the boostrap sample
												 	keep <- sample.int(nrow(Events), nrow(Events), replace = TRUE)
												 	idx <- logical( nrow(Events) )
												 	idx[keep] <- TRUE
												 	
												 	## subdivide by sex and generation
												 	Events.male.G1 <- Events[ idx & apply(outer(Events[,"Type"],c("MGP","PGP"),"=="),1,any), ]
												 	Events.female.G1 <- Events[ idx & apply(outer(Events[,"Type"],c("MGM","PGM"),"=="),1,any), ]
												 	Events.male.G2 <- Events[ idx & apply(outer(Events[,"Type"],c("P"),"=="),1,any), ]
												 	Events.female.G2 <- Events[ idx & apply(outer(Events[,"Type"],c("M"),"=="),1,any), ]
												 	
												 	## compute map length by strain, per generation
												 	G2Map.SexStrainGen <- NULL
												 	for(i in 1:8){
												 		G2Map.SexStrainGen <- rbind(G2Map.SexStrainGen,
												 																transform(padMap(ddply(Events.male.G1[Events.male.G1[,32+i],],.(Chr),integrateCrossovers)),	
												 																					Sex="Male",Gen="G1",Strain=StrainNames[i],cM=cumEvents*scalefac.A[i,1]))
												 	}
												 	for(i in 1:8){
												 		G2Map.SexStrainGen <- rbind(G2Map.SexStrainGen,
												 																transform(padMap(ddply(Events.male.G2[Events.male.G2[,32+i],],.(Chr),integrateCrossovers)),	
												 																					Sex="Male",Gen="G2",Strain=StrainNames[i],cM=cumEvents*scalefac.A[i,3]))
												 	}
												 	for(i in 1:8){
												 		G2Map.SexStrainGen <- rbind(G2Map.SexStrainGen,
												 																transform(padMap(ddply(Events.female.G1[Events.female.G1[,32+i],],
												 																											 .(Chr),integrateCrossovers)),	
												 																					Sex="Female",Gen="G1",Strain=StrainNames[i],
												 																					cM=ifelse(Chr<20,cumEvents*scalefac.A[i,2],cumEvents*scalefac.X[i,1])))
												 	}
												 	for(i in 1:8){
												 		G2Map.SexStrainGen <- rbind(G2Map.SexStrainGen,
												 																transform(padMap(ddply(Events.female.G2[Events.female.G2[,32+i],],
												 																											 .(Chr),integrateCrossovers)),	
												 																					Sex="Female",Gen="G2",Strain=StrainNames[i],
												 																					cM=ifelse(Chr<20,cumEvents*scalefac.A[i,4],cumEvents*scalefac.X[i,2])))
												 	}
												 	
												 	## return the result, with a counter indicating which replicate it is
												 	transform(G2Map.SexStrainGen, rep = j)
												 	
												 }, .progress = "text")

## total map lengths by sex, strain, gen
ddply(ddply(subset(G2Map.SexStrainGen,Chr<20),.(Chr, Sex, Gen, Strain),	
			numcolwise(max))[,c(2:4,7)],.(Strain,Sex,Gen),numcolwise(sum))

## total map lengths by sex, gen (bootstrapped for estimating variance)
g2map.summary <- ddply(ddply(subset(G2Map.bootstrap,Chr<20),
														 .(rep, Chr, Sex, Gen),	
														 colwise(max, .(cM))),
											 .(rep, Sex, Gen), colwise(sum, .(cM)))

SEX.COLORS <- setNames( c("#0080ff","#ff00ff"), levels(g2map.summary$Sex) )
ggplot(g2map.summary, aes(x = Sex, fill = Sex, y = 1.25*cM)) +
	geom_boxplot(alpha = 0.7) + 
	facet_grid(. ~ Gen) +
	scale_fill_manual(values = SEX.COLORS) + guides(fill = FALSE) +
	xlab("Generation") + ylab("autosomal map length (cM)") +
	theme_bw()

## Figure 2: scale maps by combinations of chrX and autosomal subspecific origin
Events.male.G1 <- Events[ apply(outer(Events[,"Type"],c("MGP","PGP"),"=="),1,any), ]
autosome.subspec <- apply(Events.male.G1, 1,
													function(row) {
														autosome.comb <- sort( c(row["Prdm9L.Sub"], row["Prdm9R.Sub"]) )
														paste(autosome.comb, collapse = "/")
													} )
Events.male.G1$A.strain <- autosome.subspec
scalefac.G1A.subspec <- 100/as.vector(table(Events.male.G1$A.strain))

# scaling factor for map length by chrX and autosomal subspec origin at G1 in males only
scalefac.tmp <- 100/(0.75*as.matrix( xtabs(~ X.strain + A.strain, unique(Events.male.G1[ ,c("Strain","A.strain","X.strain") ])) ))
# construct bootstrapped maps for G1 males only
g1.male.map <- ldply(1:nreps, function(i) {
	ddply(padMap(ddply(Events.male.G1[ sample.int(nrow(Events.male.G1), nrow(Events.male.G1), replace = TRUE), ],.(Chr, X.strain, A.strain), integrateCrossovers)),
				.(X.strain, A.strain),
				function(d) transform(d, Sex="Male",Gen="G1",rep=i, cM=2*cumEvents*scalefac.tmp[ unique(X.strain)[1],unique(A.strain)[1] ]) )
}, .progress = "text")

g1.male.len <- (ddply(ddply(subset(g1.male.map,Chr<20),.(rep, Chr, Sex, X.strain, A.strain),	
														numcolwise(max)), .(rep, Sex, X.strain,A.strain), colwise(sum, .(cumEvents, cM))))
g1.male.len <- transform(g1.male.len, X.strain = factor(X.strain), A.strain = factor(A.strain, levels = c("CAST/CAST","CAST/DOM","CAST/MUS","DOM/DOM","DOM/MUS","MUS/MUS")))

## use linear model to predict values for unobserved chrX/autosome combinations
g1.male.len <- transform(g1.male.len, A.CAST = as.integer(grepl("CAST", A.strain)), A.DOM = as.integer(grepl("DOM", A.strain)), A.MUS = as.integer(grepl("MUS", A.strain)))
g1.male.len <- transform(g1.male.len, A.DOM = A.DOM + ifelse(A.strain == "DOM/DOM", 1, 0), predicted = FALSE)
m0 <- lm(cM ~ X.strain + A.CAST + A.DOM + A.MUS, data = g1.male.len)
# CAST X, both CAST autosomes
predicted <- data.frame( rep = 1:100, X.strain = "CAST", A.strain = "CAST/CAST", predicted = TRUE,
												 cM = predict(m0, data.frame(X.strain = "CAST", A.CAST = 2, A.DOM = 0, A.MUS = 0)) + rnorm(100, mean(resid(m0)), sd(resid(m0))) )
# MUS X, both MUS autosomes
predicted <- rbind(predicted,
									 data.frame( rep = 1:100, X.strain = "MUS", A.strain = "MUS/MUS", predicted = TRUE,
									 						cM = predict(m0, data.frame(X.strain = "MUS", A.CAST = 0, A.DOM = 0, A.MUS = 2)) + rnorm(100, mean(resid(m0)), sd(resid(m0))) )
									 )
predicted <- rbind( g1.male.len[ ,c("rep","X.strain","A.strain","cM","predicted") ], predicted )

## Figure 2: make final plot, including both observed and predicted map lengths in G1 males
ggplot(predicted, aes(x = A.strain, y = cM)) +
	geom_boxplot(aes(fill = X.strain, alpha = !predicted, lty = predicted)) + 
	facet_grid(. ~ X.strain, scale = "free_x", drop = TRUE) +
	scale_fill_manual(values = setNames( CC.colors[ c("CAST","WSB","PWK") ], c("CAST","DOM","MUS") )) +
	guides(fill = FALSE, lty = FALSE, alpha = FALSE) +
	xlab("\nautosomal diplotype, grouped by X chromosome haplotype") +
	theme_bw() + theme( axis.text.x = element_text(angle = 90, hjust = 1) )
  
## plot cumulative sex- AND strain- AND generation-specific maps, 1 page per chromosome, to file
pdf("CumulativeBySexStrainGen.pdf", width=9.5, height=7, onefile= TRUE)		
for(i in 1:20){
	print(qplot(Pos, cM, data=subset(G2Map.SexStrainGen,Chr==i), geom="line", 
			color=Strain,main=paste("Chromosome", as.character(i), sep=" ")) +
			aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt) +
			facet_grid(Sex~Gen))
}	
dev.off()
pdf("CumulativeByStrainGenSex.pdf", width=9.5, height=7, onefile= TRUE)	
for(i in 1:20){
		print(qplot(Pos, cM, data=subset(G2Map.SexStrainGen,Chr==i), 
			geom="line", group=Sex:Gen, lty=Gen, color=Sex,
			main=paste("Chromosome", as.character(i), sep=" ")) +
			aes(color=Sex) + scale_colour_manual(value=color.sex) +
			facet_wrap(~Strain, ncol=2))
}
dev.off()

## compute recombination density
dscale = 1000
dstep = dscale/10
G2Den.SexStrainGen <- ddply(G2Map.SexStrainGen, .(Chr, Sex, Strain, Gen), 
	make.density, dstep, dscale)
head(G2Den.SexStrainGen)

## plot density to file, 1 page per chromosome
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_StrainSexGen.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexStrainGen, Chr==i), geom="line", 
		color=Strain, main=paste("Chromosome", as.character(i), sep=" ")) + 
    	aes(colour=Strain)+scale_colour_manual(value=CC.colors.alt)+
    	facet_grid(Sex~Gen))
}
dev.off()
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_StrainGenSex.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexStrainGen, Chr==i), 
		geom="line", group=Sex:Gen, lty=Gen, color=Sex,
		main=paste("Chromosome", as.character(i), sep=" ")) + 
    	aes(color=Sex)+scale_colour_manual(value=color.sex)+
    	facet_wrap(~Strain, ncol=2))
}
dev.off()



###########################################
## Compute maps by Prdm9 genotypes

Prdm9Names <- c("A","B","C","D","U")
# A is the CAST allele w/ 11 fingers
# B is B6/S129/AJ/NZO  w/ 12 fingers
# C is WSB/NOD  w/13 fingers
# D is the PWK allele with 14 fingers
Prdm9GenoNames <- c("AB","AC","AD","BA","BB","BC","BD","CA","CB","CC","CD","DA","DB","DC")

## also add a Prdm9 genotypes variable to Events
Events <- transform(Events, Prdm9Geno = paste(Events[,"Prdm9L"],Events[,"Prdm9R"],sep=''))
Events$Prdm9R <- as.character(Events$Prdm9R)
Events$Prdm9L <- as.character(Events$Prdm9L)
Events$Prdm9R[ !(Events$Prdm9R %in% Prdm9Names) ] <- "U"
Events$Prdm9L[ !(Events$Prdm9L %in% Prdm9Names) ] <- "U"
Events$Prdm9R <- factor(Events$Prdm9R, levels = Prdm9Names)
Events$Prdm9L <- factor(Events$Prdm9L, levels = Prdm9Names)
Events[which(Events$Prdm9L=="U"|Events$Prdm9R=="U"),"Prdm9Geno"] <- NaN
levels(Events$Prdm9Geno) <- Prdm9GenoNames 
table(Events$Prdm9Geno)

## total numbers of events by Prdm9 genotypes
table(Events$Prdm9L, Events$Prdm9R)

## number of meioses by Prdm9 alleles
table(Events[indx.m,]$Prdm9L, Events[indx.m,]$Prdm9R)

## avg number of events per meiosis by Prdm9 allele
round( table(Events$Prdm9L, Events$Prdm9R)[1:4,1:4] /
table(Events[indx.m,]$Prdm9L, Events[indx.m,]$Prdm9R)[1:4,1:4]  ,2)

## total events for each Prdm9 genotypes by type of meiosis
table(Events$Prdm9Geno, Events$Type)[,MeiosisNames]
table(Events$Prdm9L, Events$Type)[1:4,MeiosisNames]
table(Events$Prdm9R, Events$Type)[1:4,MeiosisNames]

## number of meioses for each Prdm9 genotypes by type of meiosis
Prdm9GenoCounts <- table(Events[indx.m,]$Prdm9Geno,Events[indx.m,]$Type)[Prdm9GenoNames,MeiosisNames]
Prdm9LCounts <- table(Events[indx.m,]$Prdm9L, Events[indx.m,]$Type)[1:4,MeiosisNames]
Prdm9RCounts <- table(Events[indx.m,]$Prdm9R, Events[indx.m,]$Type)[1:4,MeiosisNames]

## avg numbers of events
(tbl1 <- round( table(Events$Prdm9Geno, Events$Type)[,MeiosisNames]/
	table(Events[indx.m,]$Prdm9Geno, Events[indx.m,]$Type)[,MeiosisNames] ,2))
(tbl2 <- round( table(Events$Prdm9L, Events$Type)[1:4,MeiosisNames] /
	table(Events[indx.m,]$Prdm9L, Events[indx.m,]$Type)[1:4,MeiosisNames] ,2))
(tbl3 <- round( table(Events$Prdm9R, Events$Type)[1:4,MeiosisNames] /
	table(Events[indx.m,]$Prdm9R, Events[indx.m,]$Type)[1:4,MeiosisNames] ,2))
	
## compute sex ratios (female:male) by Prdm9 genotype and generation
fm.ratio <- cbind( (tbl1[,1]+tbl1[,3])/(tbl1[,2]+tbl1[,4]) , tbl1[,c(5)]/tbl1[,c(6)] )
colnames(fm.ratio)  <- c("G1", "G2")
fm.ratio

fm.Lratio <- cbind( (tbl2[,1]+tbl2[,3])/(tbl2[,2]+tbl2[,4]) , tbl2[,c(5)]/tbl2[,c(6)] )
colnames(fm.Lratio)  <- c("G1", "G2")
fm.Lratio

fm.Rratio <- cbind( (tbl3[,1]+tbl3[,3])/(tbl3[,2]+tbl3[,4]) , tbl3[,c(5)]/tbl3[,c(6)] )
colnames(fm.Rratio)  <- c("G1", "G2")
fm.Rratio

## count the frequency for any Prdm9 allele to appear in a meiosis
Prdm9Counts <- NULL
for(name in MeiosisNames) {
	Prdm9Counts <- cbind(Prdm9Counts, c(
		length(unique(subset(Events,Type==name&(Prdm9L=="A"|Prdm9R=="A"))[,"Funnel_ID"])),
		length(unique(subset(Events,Type==name&(Prdm9L=="B"|Prdm9R=="B"))[,"Funnel_ID"])),
		length(unique(subset(Events,Type==name&(Prdm9L=="C"|Prdm9R=="C"))[,"Funnel_ID"])),
		length(unique(subset(Events,Type==name&(Prdm9L=="D"|Prdm9R=="D"))[,"Funnel_ID"])),
		length(unique(subset(Events,Type==name&(Prdm9L=="U"|Prdm9R=="U"))[,"Funnel_ID"])) ))
}
rownames(Prdm9Counts) <- Prdm9Names
colnames(Prdm9Counts) <- MeiosisNames
Prdm9Counts <- Prdm9Counts[1:4,]
# (show the result)
Prdm9Counts

## count the frequency for each Prdm9 allele to NOT appear in a meiosis
Prdm9NotCounts <- NULL
for(name in MeiosisNames){
	Prdm9NotCounts <- cbind(Prdm9NotCounts, c(
		length(unique(subset(Events,
			Type==name&(Prdm9L!="A"&Prdm9R!="A"&Prdm9L!="U"&Prdm9R!="U"))[,"Funnel_ID"])),
		length(unique(subset(Events,
			Type==name&(Prdm9L!="B"&Prdm9R!="B"&Prdm9L!="U"&Prdm9R!="U"))[,"Funnel_ID"])),
		length(unique(subset(Events,
			Type==name&(Prdm9L!="C"&Prdm9R!="C"&Prdm9L!="U"&Prdm9R!="U"))[,"Funnel_ID"])),
		length(unique(subset(Events,
			Type==name&(Prdm9L!="D"&Prdm9R!="D"&Prdm9L!="U"&Prdm9R!="U"))[,"Funnel_ID"]))))
}
# this is a little sloppy wrt to handling "U", because there is 
# more information about what the allele is NOT than we have encoded here
rownames(Prdm9NotCounts) <- Prdm9Names[1:4]
colnames(Prdm9NotCounts) <- MeiosisNames
# (show the result)
Prdm9NotCounts

## look at the avg numbers of events by Meiosis and Prdm9L allele
table(Events$Prdm9L, Events$Type)[1:4,MeiosisNames]
Prdm9LCounts[1:4,]
table(Events$Prdm9L, Events$Type)[1:4,MeiosisNames]/Prdm9LCounts[1:4,]
table(Events$Prdm9R, Events$Type)[1:4,MeiosisNames]
Prdm9RCounts[1:4,]
table(Events$Prdm9R, Events$Type)[1:4,MeiosisNames]/Prdm9RCounts[1:4,]

## total number of events per Prdm9 genotype
table(Events$Prdm9L, Events$Prdm9R)[1:4,1:4]

## number of occurrences of each Prdm9 genotype
indx <- which(!duplicated(Events[,c("Funnel_ID","Type")]))
table(Events[indx,"Prdm9L"], Events[indx,"Prdm9R"])[1:4,1:4]
table(Events$Prdm9L, Events$Prdm9R)[1:4,1:4]/table(Events[indx,"Prdm9L"], Events[indx,"Prdm9R"])[1:4,1:4]

tmp <- matrix(c(3/4,3/4,3/4,3/4,2,2),ncol=1)
ScaleFac.A <- 100/Prdm9Counts%*%tmp

tmp <- matrix(c(3/4,0,1,0,2,0),ncol=1)
ScaleFac.X <- 100/Prdm9Counts%*%tmp

## compute maps per Prdm9 allele
G2Map.Prdm9 <- rbind(
  transform(padMap(ddply(subset(Events,Prdm9L=="A"|Prdm9R=="A"),.(Chr),integrateCrossovers)),
  		Allele="Prdm9_A",
  		cM=ifelse(Chr<20, cumEvents*ScaleFac.A[1], cumEvents*ScaleFac.X[1])),  
   transform(padMap(ddply(subset(Events,Prdm9L=="B"|Prdm9R=="B"),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_B",
   		cM=ifelse(Chr<20, cumEvents*ScaleFac.A[2], cumEvents*ScaleFac.X[2])),  
   transform(padMap(ddply(subset(Events,Prdm9L=="C"|Prdm9R=="C"),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_C",
   		cM=ifelse(Chr<20, cumEvents*ScaleFac.A[3], cumEvents*ScaleFac.X[3])),  
   transform(padMap(ddply(subset(Events,Prdm9L=="D"|Prdm9R=="D"),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_D",
   		cM=ifelse(Chr<20, cumEvents*ScaleFac.A[4], cumEvents*ScaleFac.X[4])))  
head(G2Map.Prdm9)

## check scaling
ddply(ddply(G2Map.Prdm9,.(Chr, Allele),numcolwise(max)),
	.(Chr),numcolwise(mean))$cM / ddply(G2Map.avg,.(Chr),numcolwise(max))$cM
	
## plot cumulative maps by Prdm9 allele, 1 chromosome per page, to file
pdf("CumulativeByPrdm9.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
		print(qplot(Pos, cM, data=subset(G2Map.Prdm9,Chr==i), geom="line", 
	color=Allele,main=paste("Chromosome", as.character(i), sep=" ")))
	}
dev.off()

## compute recombination density by Prdm9 allele
dscale = 1000
dstep = dscale/10
G2Den.Prdm9 <- ddply(G2Map.Prdm9, .(Chr, Allele), make.density, dstep, dscale)
head(G2Den.Prdm9)

## plot it to file, 1 chromosome per page
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_Prdm9.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.Prdm9, Chr==i), geom="line", color=Allele,
		main=paste("Chromosome", as.character(i), sep=" ")) +
	    facet_grid(facets=Allele ~.) )
}
dev.off()



####################################################
## Compute sex- and Prdm9 allele-specific maps in G2F1 population

Events.male <- Events[apply(outer(Events[,"Type"],c("MGP","PGP","P"),"=="),1,any),]
Events.female <- Events[apply(outer(Events[,"Type"],c("MGM","PGM","M"),"=="),1,any),]

## prepare scaling factors
ScaleFac.A <- 100/cbind((Prdm9Counts[,1]+Prdm9Counts[,3])*3/4 + Prdm9Counts[,5]*2,
	(Prdm9Counts[,2]+Prdm9Counts[,4])*3/4 + Prdm9Counts[,6]*2)
colnames(ScaleFac.A) <- c("Female","Male")
ScaleFac.X <- 100/cbind(Prdm9Counts[,1]*3/4+Prdm9Counts[,3] + Prdm9Counts[,5]*2,c(0,0,0,0))
colnames(ScaleFac.X) <- c("Female","Male")

## compute integrated maps
G2Map.SexPrdm9 <- rbind(
  transform(padMap(ddply(subset(Events.male,Prdm9L=="A"|Prdm9R=="A"),
  		.(Chr),integrateCrossovers)), Allele="Prdm9_A",Sex="Male",
  		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["A","Male"],cumEvents*ScaleFac.X["A","Male"])),  
   transform(padMap(ddply(subset(Events.male,Prdm9L=="B"|Prdm9R=="B"),
   		.(Chr),integrateCrossovers)), Allele="Prdm9_B",Sex="Male",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["B","Male"],cumEvents*ScaleFac.X["B","Male"])),  
   transform(padMap(ddply(subset(Events.male,Prdm9L=="C"|Prdm9R=="C"),
   		.(Chr),integrateCrossovers)), Allele="Prdm9_C",Sex="Male",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["C","Male"],cumEvents*ScaleFac.X["C","Male"])),  
   transform(padMap(ddply(subset(Events.male,Prdm9L=="D"|Prdm9R=="D"),
   		.(Chr),integrateCrossovers)), Allele="Prdm9_D",Sex="Male",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["D","Male"],cumEvents*ScaleFac.X["D","Male"])),
   transform(padMap(ddply(subset(Events.female,Prdm9L=="A"|Prdm9R=="A"),
   		.(Chr),integrateCrossovers)), Allele="Prdm9_A",Sex="Female",
  		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["A","Female"],cumEvents*ScaleFac.X["A","Female"])),  
   transform(padMap(ddply(subset(Events.female,Prdm9L=="B"|Prdm9R=="B"),
   		.(Chr),integrateCrossovers)), Allele="Prdm9_B",Sex="Female",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["B","Female"],cumEvents*ScaleFac.X["B","Female"])),  
   transform(padMap(ddply(subset(Events.female,Prdm9L=="C"|Prdm9R=="C"),
   		.(Chr),integrateCrossovers)), Allele="Prdm9_C",Sex="Female",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["C","Female"],cumEvents*ScaleFac.X["C","Female"])),  
   transform(padMap(ddply(subset(Events.female,Prdm9L=="D"|Prdm9R=="D"),
   		.(Chr),integrateCrossovers)), Allele="Prdm9_D",Sex="Female",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["D","Female"],cumEvents*ScaleFac.X["D","Female"])))  
# (peek to check result)
head(G2Map.SexPrdm9)

## check scaling
ddply(ddply(G2Map.SexPrdm9,.(Chr, Allele, Sex),numcolwise(max)),
	.(Chr),numcolwise(mean))$cM / ddply(G2Map.avg,.(Chr),numcolwise(max))$cM

## update chromosome lengths
ChrLen <- transform(ChrLen,
	Prdm9_A_Male_cM = c(ddply(subset(G2Map.SexPrdm9,Allele=="Prdm9_A"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_B_Male_cM = c(ddply(subset(G2Map.SexPrdm9,Allele=="Prdm9_B"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_C_Male_cM = c(ddply(subset(G2Map.SexPrdm9,Allele=="Prdm9_C"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_D_Male_cM = c(ddply(subset(G2Map.SexPrdm9,Allele=="Prdm9_D"&Sex=="Male"),.(Chr),	
		numcolwise(max))[,"cM"],NA),
	Prdm9_A_Female_cM = ddply(subset(G2Map.SexPrdm9,Allele=="Prdm9_A"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_B_Female_cM = ddply(subset(G2Map.SexPrdm9,Allele=="Prdm9_B"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_C_Female_cM = ddply(subset(G2Map.SexPrdm9,Allele=="Prdm9_C"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_D_Female_cM = ddply(subset(G2Map.SexPrdm9,Allele=="Prdm9_D"&Sex=="Female"),.(Chr),	
		numcolwise(max))[,"cM"]
)
ChrLen 

## plot cumulative sex- and Prdm9-allele specific maps, 1 page per chromosme, to file
pdf("CumulativePrdm9Sex.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
		print(qplot(Pos, cM, data=subset(G2Map.SexPrdm9,Chr==i), geom="line", lty=Sex,
	color=Allele,main=paste("Chromosome", as.character(i), sep=" ")) +
	facet_grid(Sex~.)	)
}
dev.off()
pdf("CumulativeSexPrdm9.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
		print(qplot(Pos, cM, data=subset(G2Map.SexPrdm9,Chr==i), geom="line", lty=Sex,
	color=Allele,main=paste("Chromosome", as.character(i), sep=" ")) +
	facet_grid(Allele~.)	)
	}
dev.off()

## compute recombination density by sex and Prdm9 allele
dscale = 1000
dstep = dscale/10
G2Den.SexPrdm9 <- ddply(G2Map.SexPrdm9, .(Chr, Sex, Allele), make.density, dstep, dscale)
head(G2Den.SexPrdm9)

## plot density to file, 1 page per chromosome
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_Prdm9Sex.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexPrdm9, Chr==i), geom="line", linetype=Sex,
		colour=Allele, main=paste("Chromosome", as.character(i), sep=" ")) + 
    	facet_wrap(facets=~Allele, nrow=4) )
    }
dev.off()
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_SexPrdm9.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexPrdm9, Chr==i), geom="line",
		colour=Allele, main=paste("Chromosome", as.character(i), sep=" ")) + 
    	facet_grid(Sex~.))
    }
dev.off()
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_SexPrdm92.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexPrdm9, Chr==i), geom="line",
		colour=Allele, main=paste("Chromosome", as.character(i), sep=" ")) + 
    	facet_grid(Allele~Sex))
    }
dev.off()



###########################################
## Compute maps by Prdm9 genotype, sex and generation in G2F1 population

Events <- transform(Events, Prdm9_A=Prdm9L=="A"|Prdm9R=="A", 
	Prdm9_B=Prdm9L=="B"|Prdm9R=="B", 
	Prdm9_C=Prdm9L=="C"|Prdm9R=="C", 
	Prdm9_D=Prdm9L=="D"|Prdm9R=="D" )
	
Events.male.G1 <- Events[apply(outer(Events[,"Type"],c("MGP","PGP"),"=="),1,any),]
Events.female.G1 <- Events[apply(outer(Events[,"Type"],c("MGM","PGM"),"=="),1,any),]
Events.male.G2 <- Events[apply(outer(Events[,"Type"],c("P"),"=="),1,any),]
Events.female.G2 <- Events[apply(outer(Events[,"Type"],c("M"),"=="),1,any),]
	
## prepare scaling factors
ScaleFac.A <- 100/cbind(
	(Prdm9Counts[,1]+Prdm9Counts[,3])*3/4,
	(Prdm9Counts[,2]+Prdm9Counts[,4])*3/4,
	Prdm9Counts[,5]*2, Prdm9Counts[,6]*2)
colnames(ScaleFac.A) <- c("G1_Female","G1_Male","G2_Female","G2_Male")
ScaleFac.X <- 100/cbind(
	Prdm9Counts[,1]*3/4+Prdm9Counts[,3],
	0,Prdm9Counts[,5]*2, 0)
colnames(ScaleFac.X) <- c("G1_Female","G1_Male","G2_Female","G2_Male")

## compute integrated maps
G2Map.SexGenPrdm9 <- rbind(
  transform(padMap(ddply(subset(Events.male.G1,Prdm9_A),.(Chr),integrateCrossovers)),
  		Allele="Prdm9_A",Sex="Male",Gen="G1",
  		cM=cumEvents*ScaleFac.A["A","G1_Male"]),  
   transform(padMap(ddply(subset(Events.male.G1,Prdm9_B),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_B",Sex="Male",Gen="G1",
   		cM=cumEvents*ScaleFac.A["B","G1_Male"]),  
   transform(padMap(ddply(subset(Events.male.G1,Prdm9_C),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_C",Sex="Male",Gen="G1",
   		cM=cumEvents*ScaleFac.A["C","G1_Male"]),  
   transform(padMap(ddply(subset(Events.male.G1,Prdm9_D),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_D",Sex="Male",Gen="G1",
   		cM=cumEvents*ScaleFac.A["D","G1_Male"]),
   transform(padMap(ddply(subset(Events.female.G1,Prdm9_A),.(Chr),integrateCrossovers)),
  		Allele="Prdm9_A",Sex="Female",Gen="G1",
  		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["A","G1_Female"],cumEvents*ScaleFac.X["A","G1_Female"])),  
   transform(padMap(ddply(subset(Events.female.G1,Prdm9_B),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_B",Sex="Female",Gen="G1",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["B","G1_Female"],cumEvents*ScaleFac.X["B","G1_Female"])),  
   transform(padMap(ddply(subset(Events.female.G1,Prdm9_C),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_C",Sex="Female",Gen="G1",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["C","G1_Female"],cumEvents*ScaleFac.X["C","G1_Female"])),  
   transform(padMap(ddply(subset(Events.female.G1,Prdm9_D),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_D",Sex="Female",Gen="G1",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["D","G1_Female"],cumEvents*ScaleFac.X["D","G1_Female"])),
   transform(padMap(ddply(subset(Events.male.G2,Prdm9_A),.(Chr),integrateCrossovers)),
  		Allele="Prdm9_A",Sex="Male",Gen="G2",
  		cM=cumEvents*ScaleFac.A["A","G2_Male"]),  
   transform(padMap(ddply(subset(Events.male.G2,Prdm9_B),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_B",Sex="Male",Gen="G2",
   		cM=cumEvents*ScaleFac.A["B","G2_Male"]),  
   transform(padMap(ddply(subset(Events.male.G2,Prdm9_C),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_C",Sex="Male",Gen="G2",
   		cM=cumEvents*ScaleFac.A["C","G2_Male"]),  
   transform(padMap(ddply(subset(Events.male.G2,Prdm9_D),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_D",Sex="Male",Gen="G2",
   		cM=cumEvents*ScaleFac.A["D","G2_Male"]),
   transform(padMap(ddply(subset(Events.female.G2,Prdm9_A),.(Chr),integrateCrossovers)),
  		Allele="Prdm9_A",Sex="Female",Gen="G2",
  		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["A","G2_Female"],cumEvents*ScaleFac.X["A","G2_Female"])),  
   transform(padMap(ddply(subset(Events.female.G2,Prdm9_B),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_B",Sex="Female",Gen="G2",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["B","G2_Female"],cumEvents*ScaleFac.X["B","G2_Female"])),  
   transform(padMap(ddply(subset(Events.female.G2,Prdm9_C),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_C",Sex="Female",Gen="G2",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["C","G2_Female"],cumEvents*ScaleFac.X["C","G2_Female"])),  
   transform(padMap(ddply(subset(Events.female.G2,Prdm9_D),.(Chr),integrateCrossovers)),
   		Allele="Prdm9_D",Sex="Female",Gen="G2",
   		cM=ifelse(Chr<20,cumEvents*ScaleFac.A["D","G2_Female"],cumEvents*ScaleFac.X["D","G2_Female"])))  
# (peek to check result)
head(G2Map.SexGenPrdm9)

## check scaling
ddply(ddply(G2Map.SexGenPrdm9,.(Chr, Allele, Sex, Gen),numcolwise(max)),
	.(Chr),numcolwise(mean))$cM / ddply(G2Map.avg,.(Chr),numcolwise(max))$cM\
ddply(ddply(subset(G2Map.SexGenPrdm9,Gen=="G1"),.(Chr,Sex,Allele),numcolwise(max)),
	.(Chr),numcolwise(mean))$cM / ddply(G2Map.avg,.(Chr),numcolwise(max))$cM 
ddply(ddply(subset(G2Map.SexGenPrdm9,Gen=="G2"),.(Chr,Sex,Allele),numcolwise(max)),
	.(Chr),numcolwise(mean))$cM / ddply(G2Map.avg,.(Chr),numcolwise(max))$cM  

## update chromosome lengths
ChrLen <- transform(ChrLen,
	Prdm9_A_Male_Gen1_cM = c(ddply(subset(G2Map.SexGenPrdm9,	
		Allele=="Prdm9_A"&Sex=="Male"&Gen=="G1"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_B_Male_Gen1_cM = c(ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_B"&Sex=="Male"&Gen=="G1"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_C_Male_Gen1_cM = c(ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_C"&Sex=="Male"&Gen=="G1"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_D_Male_Gen1_cM = c(ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_D"&Sex=="Male"&Gen=="G1"),.(Chr),	
		numcolwise(max))[,"cM"],NA),
	Prdm9_A_Female_Gen1_cM = ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_A"&Sex=="Female"&Gen=="G1"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_B_Female_Gen1_cM = ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_B"&Sex=="Female"&Gen=="G1"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_C_Female_Gen1_cM = ddply(subset(G2Map.SexGenPrdm9,	
		Allele=="Prdm9_C"&Sex=="Female"&Gen=="G1"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_D_Female_Gen1_cM = ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_D"&Sex=="Female"&Gen=="G1"),.(Chr),	
		numcolwise(max))[,"cM"],
	Prdm9_A_Male_Gen2_cM = c(ddply(subset(G2Map.SexGenPrdm9,	
		Allele=="Prdm9_A"&Sex=="Male"&Gen=="G2"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_B_Male_Gen2_cM = c(ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_B"&Sex=="Male"&Gen=="G2"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_C_Male_Gen2_cM = c(ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_C"&Sex=="Male"&Gen=="G2"),.(Chr),	
		numcolwise(max))[,"cM"],NA),	
	Prdm9_D_Male_Gen2_cM = c(ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_D"&Sex=="Male"&Gen=="G2"),.(Chr),	
		numcolwise(max))[,"cM"],NA),
	Prdm9_A_Female_Gen2_cM = ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_A"&Sex=="Female"&Gen=="G2"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_B_Female_Gen2_cM = ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_B"&Sex=="Female"&Gen=="G2"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_C_Female_Gen2_cM = ddply(subset(G2Map.SexGenPrdm9,	
		Allele=="Prdm9_C"&Sex=="Female"&Gen=="G2"),.(Chr),	
		numcolwise(max))[,"cM"],	
	Prdm9_D_Female_Gen2_cM = ddply(subset(G2Map.SexGenPrdm9,
		Allele=="Prdm9_D"&Sex=="Female"&Gen=="G2"),.(Chr),	
		numcolwise(max))[,"cM"]
)
ChrLen 

## plot cumulative maps by sex, generation and Prdm9 allele, 1 page per chromosome, to file
pdf("CumulativeBySexGenPrdm9.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
		print(qplot(Pos, cM, data=subset(G2Map.SexGenPrdm9,Chr==i), geom="line", lty=Sex,
	color=Allele,main=paste("Chromosome", as.character(i), sep=" ")) +
	facet_grid(Sex~Gen)	)
}
dev.off()
pdf("CumulativeByPrdm9SexGen.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
		print(qplot(Pos, cM, data=subset(G2Map.SexGenPrdm9,Chr==i), geom="line", lty=Sex,
	color=Allele,main=paste("Chromosome", as.character(i), sep=" ")) +
	facet_grid(Allele~Gen)  )
}	
dev.off()
pdf("CumulativeByGenPrdm9Sex.pdf", width=9.5, height=7, onefile= TRUE)
for(i in 1:20){
		print(qplot(Pos, cM, data=subset(G2Map.SexGenPrdm9,Chr==i), geom="line", lty=Gen,
	color=Allele,main=paste("Chromosome", as.character(i), sep=" ")) +
	facet_grid(Allele~Sex)	)
}
dev.off()

## compute recombination density by sex, generation and Prdm9 allele
dscale = 1000
dstep = dscale/10
G2Den.SexGenPrdm9 <- ddply(G2Map.SexGenPrdm9, .(Chr, Sex, Allele, Gen), 
	make.density, dstep, dscale)
head(G2Den.SexGenPrdm9)

## plot density to file, 1 page per chromosome
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_SexGenSex.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexGenPrdm9, Chr==i), geom="line", 
		 color=Allele, main=paste("Chromosome", as.character(i), sep=" ")) + 
    	 facet_grid(Sex~Gen))
}
dev.off()
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_Prdm9SexGen.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexGenPrdm9, Chr==i), 
		geom="line", color=Sex,
		main=paste("Chromosome", as.character(i), sep=" ")) + 
    	aes(color=Sex)+scale_colour_manual(value=color.sex)+
    	facet_grid(Allele~Gen))
}
dev.off()
pdf(paste("Density_G2F1_",as.character(dscale/1000),"Mb_Prdm9GenSex.pdf",sep=""))
for(i in 1:20){
	print(qplot(Pos,cM,data=subset(G2Den.SexGenPrdm9, Chr==i), 
		geom="line", color=Gen,
		main=paste("Chromosome", as.character(i), sep=" ")) + 
    	facet_grid(Allele~Sex))
    }
dev.off()



##########################################
## write chromosome lengths to file
write.csv(ChrLen,"ChrLen.csv")
##########################################