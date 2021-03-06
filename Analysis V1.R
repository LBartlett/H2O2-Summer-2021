##

# Analysis for Manuscript:

# Contact: lewis.bartlett@uga.edu    https://scholar.google.com/citations?user=PV5ca5UAAAAJ&hl=en

## 

library(afex)

library(emmeans)

## Read in Data Sets

BMort <- read.csv(file = 'BMort.csv',
                  header = TRUE)

SMort <- read.csv(file = 'SMort.csv',
                  header = TRUE)

S50 <- read.csv(file = 'S50.csv',
                header = TRUE)

G50 <- read.csv(file = 'G50.csv',
                header = TRUE)

SLong <- read.csv(file = 'SLong.csv',
                  header = TRUE)

GLong <- read.csv(file = 'GLong.csv',
                  header = TRUE)

###
# Mortality analysis first

# Lower concentrations

SMort$DA <- cbind(SMort$Died, SMort$Alive)

SMort$Cage <- (paste0((1:NROW(SMort)),SMort$Colony))

MortMod <- mixed(DA ~ Concentration + (1|Colony) + (1|Colony:Cage),
                 family = 'binomial',
                 data = SMort,
                 method = 'LRT')

nice(MortMod)

summary(MortMod)

# No strong evidence for any increasing mortality with changing concentration at these lower concentrations


# Higher concentrations 

BMort$DA <- cbind(BMort$Dead, BMort$Alive)

MortModB <- mixed(DA ~ Dose + Starved + (1|Colony) + (1|Colony:Pot),
                  family = 'binomial',
                  data = BMort,
                  method = 'LRT')

nice(MortModB)

summary(MortModB)

# Strong evidence of increasing mortality with increasing concentration at these higher concentrations


## Choice Analysis

# 50ug/l

t.test(G50$RV, mu = 0, alternative = 'two.sided')
t.test(S50$RV, mu = 0, alternative = 'two.sided')

C50Full <- rbind(G50,S50)

CMod50 <- mixed(RV ~ Sugar + (1|Sugar:UCol),
                data = C50Full,
                test_intercept = TRUE)

nice(CMod50)

summary(CMod50)

emmeans(CMod50, specs = 'Sugar')

# Evidence of *avoiding* H202 at this concentration, with effect size  stronger in glucose, but not significantly.


# Concentration-choice

SModLong <- mixed(RV ~ Concentration + (1|UCol),
                  data = SLong,
                  test_intercept = TRUE)
nice(SModLong)

# No good evidence of H202 avoidance in sucrose here

GModLong <- mixed(RV ~ Concentration + (1|UCol),
                  data = GLong,
                  test_intercept = TRUE)
nice(GModLong)

summary(GModLong)

# Mixed evidence of H2O2 avoidance in glucose

CLongFull <- rbind(GLong,SLong)

CModLong <- mixed(RV ~ Concentration*Sugar + (1|Sugar:UCol),
                  data = CLongFull,
                  test_intercept = TRUE)

nice(CModLong)

summary(CModLong)

# Avoidance present, but no difference betwene sugars and a lessening of avoidance at higher concentrations?



#### Graphs

# Main mortality

PMod <- glmer(DA ~ Dose + Starved + (1|Colony),
              family = 'binomial',
              data = BMort)

par(mar = c(5,5,5,2))

BMort$PD <- BMort$Dead/BMort$Total

cols <- c('red2', 'blue2', 'purple')

plot(x=(0:10), 
     y = (0:10)*0.1, 
     type='n', 
     xlab = expression(paste('Dose (% H'[2],'O'[2],')',sep='')), 
     ylab = 'Proportion Dead', 
     main = 'Bee Mortality',
     cex.main = 1.75, cex.lab = 1.635, cex.axis = 1.635)


for(OriginHive in unique(BMort$Colony)){
  
  GV1 <- seq(0,10, length.out = 500)
  
  GV2 <- predict(PMod, 
                 newdata = data.frame(Dose = GV1, 
                                      Starved = FALSE, 
                                      Colony = rep(OriginHive, times = NROW(GV1))
                 ),
                 type = 'response')
  
  points(x = GV1, y = GV2, type = 'l', lwd = 4 ,col = cols[which(unique(BMort$Colony)==OriginHive)], pch = 20, cex=1.75)
  
  points(x = jitter(BMort$Dose[which(BMort$Colony == OriginHive & BMort$Starved == F)], factor = 0.5), 
         y = BMort$PD[which(BMort$Colony == OriginHive & BMort$Starved == F)], 
         type = 'p', lwd = 2 ,
         col = cols[which(unique(BMort$Colony)==OriginHive)], 
         pch = 20, cex=2.5)
  
}


## Supplemental mortality

par(mar = c(5,5,5,2))

SMort$PD <- SMort$Died/SMort$Total

cols <- c('pink3', 'orange1', 'lightblue')

plot(x=log10(10^(((0:8)/2))), 
     y = (0:8)*0.125, 
     type='n', 
     xlab = expression(paste('Dose (log'[10],'+1 H'[2],'O'[2],')',sep='')), 
     ylab = 'Proportion Dead', 
     main = 'Bee Mortality',
     cex.main = 1.75, cex.lab = 1.635, cex.axis = 1.635)


for(OriginHive in unique(SMort$Colony)){
  
  points(x = jitter(log10(SMort$Concentration[which(SMort$Colony == OriginHive)]+1), factor = 0.15), 
         y = SMort$PD[which(SMort$Colony == OriginHive)], 
         type = 'p', lwd = 2 ,
         col = cols[which(unique(SMort$Colony)==OriginHive)], 
         pch = 20, cex=2.5)
  
}


# Choice graphs

# long choice

SLong$LConc <- log10(SLong$Concentration)
GLong$LConc <- log10(GLong$Concentration)

par(mfrow = c(2,1),
    mar = c(5,5,4,1))

BP.SLONG <- boxplot(SLong$RV ~ SLong$LConc, 
                    main = 'Sucrose', 
                    xlab = expression(paste('log'[10],' H'[2],'O'[2],' Concentration (',mu,'g ',ml^-1,')]',sep='')),
                    ylab = 'Feeding bias (ml)', 
                    border = 'transparent', cex.axis = 1.3, cex.lab = 1.5, outline = FALSE, lwd = 1.2,
                    boxlty = 0, whisklty = 0, staplelty = 0, boxwex = 0)

stripchart(SLong$RV ~ SLong$LConc,
           #col = sapply(PlotCols, Transpa, percent = 65),
           vertical = T, add = T, pch = 4, cex = 1.5, 
           method = 'jitter', lwd = 2)

PlotAg <- aggregate(SLong$RV ~ SLong$LConc, FUN = mean)
stripchart(PlotAg[,2] ~ PlotAg[,1],
           #col = sapply(PlotCols, Transpa, percent = 0),
           vertical = T, add = T, pch = 95, cex = 4, lwd = 2)

se <- function(x) sqrt(var(x)/length(x))
PAGSE <-  aggregate(SLong$RV ~ SLong$LConc, FUN = sd)
PAGSE$Lower <- PlotAg$`SLong$RV` - PAGSE$`SLong$RV`
PAGSE$Upper <- PlotAg$`SLong$RV` + PAGSE$`SLong$RV`

for(J in 1:6){
  
  segments(x0 = J, x1 = J, y0 = PAGSE$Lower[J], y1 = PAGSE$Upper[J],
           #col = PlotCols[J],
           lwd = 3)
  
}

abline(a = 0, b = 0, lty = 2)


BP.GLONG <- boxplot(GLong$RV ~ GLong$LConc, 
                    main = 'Glucose', 
                    xlab = expression(paste('log'[10],'[H'[2],'O'[2],' Concentration (',mu,'g ',ml^-1,')]',sep='')),
                    ylab = 'Feeding bias (ml)', 
                    border = 'transparent', cex.axis = 1.3, cex.lab = 1.5, outline = FALSE, lwd = 1.2,
                    boxlty = 0, whisklty = 0, staplelty = 0, boxwex = 0)

stripchart(GLong$RV ~ GLong$LConc,
           #col = sapply(PlotCols, Transpa, percent = 65),
           vertical = T, add = T, pch = 4, cex = 1.5, 
           method = 'jitter', lwd = 2)

PlotAg <- aggregate(GLong$RV ~ GLong$LConc, FUN = mean)
stripchart(PlotAg[,2] ~ PlotAg[,1],
           #col = sapply(PlotCols, Transpa, percent = 0),
           vertical = T, add = T, pch = 95, cex = 4, lwd = 2)

se <- function(x) sqrt(var(x)/length(x))
PAGSE <-  aggregate(GLong$RV ~ GLong$LConc, FUN = sd)
PAGSE$Lower <- PlotAg$`GLong$RV` - PAGSE$`GLong$RV`
PAGSE$Upper <- PlotAg$`GLong$RV` + PAGSE$`GLong$RV`

for(J in 1:6){
  
  segments(x0 = J, x1 = J, y0 = PAGSE$Lower[J], y1 = PAGSE$Upper[J],
           #col = PlotCols[J],
           lwd = 3)
  
}
abline(a = 0, b = 0, lty = 2)


# 50

par(mfrow = c(1,1))

BP.50 <- boxplot(C50Full$RV ~ C50Full$Sugar, 
                    main = NA, 
                    xlab = expression(paste('Sugar Solution ( 50',mu,'g ',ml^-1,'  H'[2],'O'[2],')',sep='')),
                    ylab = 'Feeding bias (ml)', 
                    border = 'transparent', cex.axis = 1.3, cex.lab = 1.5, outline = FALSE, lwd = 1.2,
                    boxlty = 0, whisklty = 0, staplelty = 0, boxwex = 0)

stripchart(C50Full$RV ~ C50Full$Sugar,
           #col = sapply(PlotCols, Transpa, percent = 65),
           vertical = T, add = T, pch = 4, cex = 1.5, 
           method = 'jitter', lwd = 2)

PlotAg <- aggregate(C50Full$RV ~ C50Full$Sugar, FUN = mean)
stripchart(PlotAg[,2] ~ PlotAg[,1],
           #col = sapply(PlotCols, Transpa, percent = 0),
           vertical = T, add = T, pch = 95, cex = 4, lwd = 2)

se <- function(x) sqrt(var(x)/length(x))
PAGSE <-  aggregate(C50Full$RV ~ C50Full$Sugar, FUN = sd)
PAGSE$Lower <- PlotAg$`C50Full$RV` - PAGSE$`C50Full$RV`
PAGSE$Upper <- PlotAg$`C50Full$RV` + PAGSE$`C50Full$RV`

for(J in 1:2){
  
  segments(x0 = J, x1 = J, y0 = PAGSE$Lower[J], y1 = PAGSE$Upper[J],
           #col = PlotCols[J],
           lwd = 3)
  
}
abline(a = 0, b = 0, lty = 2)



# Review Recommendation (constructive, appreciated)
# Check if consumption in no-choice is determined by concentrations
# Only really meaningful with the SMort data where mortality didn't change across concentrations

plot(SMort$Consumed[which(!SMort$Starved)]~log(SMort$Concentration[which(!SMort$Starved)]+1))

#think it's obvious but check anyway

hist(SMort$Consumed[which(!SMort$Starved)])

#declining ti pursue any stats as it's an odd measure (not ad libitum)








