

# R script to perform species distribution modelling using glm and prepare prediction maps (species distribution maps)
# Ugyen Penjor, 2019 (ugyenpenjor.bt@gmail.com)
# Link to the paper: https://doi.org/10.1007/s10980-021-01225-7
# Penjor, U., Kaszta, Z., Macdonald, D. W., & Cushman, S. A. (2021) 
# Prioritizing areas for conservation outside the existing protected area network in Bhutan: the use of multi-species, multi-scale habitat suitability models

# Set working directory
setwd("your working directory")
getwd()

# Load data
det.data <- read.csv("glm_species_distribution_model_covariates.csv", header=T)
head(det.data)

# Load presence/absence data and other input files for prediction
# Note that this procedure will only produce apparent distribution map because we are not accounting for detection probability.
load("glm_species_distribution_model_data.RData") # heavy object!
ls() # check the objects inside this packaged data

# Species presence/absence points (but we don't have temporal replicates)
head(y.det)

#############################################################
########## Check for spatial autocorrelation (SAC) ##########
#############################################################

# fit a logistic regression for naive occupancy
( mod <- glm(y.det ~ scale(AGRI_GYRATE_1KM) + scale(BARE_GYRATE_2KM) + scale(BLF_PD_32km) + scale(BUILTUP_GYR_4km) +scale(LC_GYRATE_1km) + 
               scale(LC_CWED_4km) + scale(TC_CL2_GYRATE_8km) + scale(rivFM_2KM) + scale(roaFM_1KM) + scale(sloposFM_4KM) + scale(effort), 
             data= det.data, family=binomial(link="logit")) )

# Load the package to use the function
library(ncf)

# Spline correlogram for naive occupancy 
Correlog.mod <- spline.correlog(
  x=det.data[,2],
  y=det.data[,3],
  z=residuals(mod, type="pearson"), xmax=F,
  latlon=T) # if xmax=F takes time to run, if xmax=25, short run

plot(Correlog.mod, ylim=c(-0.35, 0.35), ylab="Correlation (Moran's I)", xlab="Distance (m)", main="Spatial correlogram for Asiatic wild dog")

summary(Correlog.mod)
print.default(Correlog.mod)

# Spatial autocorrelation (SAC) test 
library(spdep)

# We need species presence points and coordinates
z <- det.data$dhole
x <- det.data$X
y <- det.data$Y
xy <- cbind(x,y)

# Nearest neighbhour analysis
pts_kn1 <- knn2nb(knearneigh(xy, k=1))

# Distance between observation points (in this case camera stations)
dist <- unlist(nbdists(pts_kn1, xy))

( d_mn <- mean(dist) )  # Mean distance between points
( d_sd <- sd(dist) )    # standard deviation

# We're only concerned with autocorrelation WITHIN sites, not among. Maybe max of within-site mean dist.?
d_mn <- 2910.584 # mean distance between camera traps in metres

pts_kd_mn <- dnearneigh(xy, d1=0, d2=d_mn)  #d1 = starting distance; d2 is max distance to consider for spatial autocorrelation
pts_kd_mn_w <- nb2listw(pts_kd_mn, zero.policy=TRUE)

moran.test(z, pts_kd_mn_w, zero.policy=TRUE)

# Moran's I test to check SAC
test <- glm(y.det ~ scale(AGRI_GYRATE_1KM) + scale(BARE_GYRATE_2KM) + scale(BLF_PD_32km) + scale(BUILTUP_GYR_4km) +scale(LC_GYRATE_1km) + 
              scale(LC_CWED_4km) + scale(TC_CL2_GYRATE_8km) + scale(rivFM_2KM) + scale(roaFM_1KM) + scale(sloposFM_4KM) + scale(effort), 
            data= det.data, family=binomial(link="logit"), na.action="na.fail") 

summary(test)

lm.morantest(test, pts_kd_mn_w, zero.policy=T) 

## lm doesn't work for the glmer mixed effects model
#vars <- cbind(file$CF_CorrLength_Scale500.s,
#              file$NF_CorrLength_Scale4000.s, 
#              file$OF_CorrLength_Scale10000.s, 
#              file$CTI_FM_Scale8000.s, 
#              file$CTI_SD_Scale8000.s,
#              file$CTI_SD_Scale8000.sq.s,
#              file$TreesReclass_ED_Scale1000.s)

#moran.test(vars, pts_kd_mn_w, zero.policy=TRUE)
#length(vars)
#length(pts_kd_mn_w)

# If p-value is significant, run the following code, otherwise can skip

# Prepare an autocovariate to correct for predictions among the nearest neighbhours
( aut <- autocov_dist(z, xy, nbs=d_mn, type="inverse", zero.policy=T) )

# If aut is significant, add the aut term to the model below

mod.aut <- glm(y.det ~scale(AGRI_GYRATE_1KM) + scale(BARE_GYRATE_2KM) + scale(BLF_PD_32km) + scale(BUILTUP_GYR_4km) +scale(LC_GYRATE_1km) + 
                 scale(LC_CWED_4km) + scale(TC_CL2_GYRATE_8km) + scale(rivFM_2KM) + scale(roaFM_1KM) + scale(sloposFM_4KM) + scale(effort) +
                 scale(effort) + aut, # autocovariate added here
               data= det.data, family=binomial(link="logit"), na.action="na.fail") 

lm.morantest(mod.aut, pts_kd_mn_w, zero.policy=T)

# New moran test with aut term included (it should correct for the SAC, if not, we need to test with different distance value and re-run it again)

summary(mod.aut)


# Multivariate modelling - here we did not find any significant affect of SAC, so we continue with the model that doesn't account for SAC
# Model without autocovariate
(glob.mod <- glm(y.det ~ scale(AGRI_GYRATE_1KM) + scale(BARE_GYRATE_2KM) + scale(BLF_PD_32km) + scale(BUILTUP_GYR_4km) +scale(LC_GYRATE_1km) + 
                   scale(LC_CWED_4km) + scale(TC_CL2_GYRATE_8km) + scale(rivFM_2KM) + scale(roaFM_1KM) + scale(sloposFM_4KM) + scale(effort),
                 data= det.data, family=binomial(link="logit"))) 

library(MuMIn)
# All possible combination of predictors
options(na.action=na.fail)  # to handle missing values

(fm_dredge <- dredge(glob.mod, rank="AICc", trace=TRUE))

# Sum of Aikake weights (Summed Model Weights (SMW))
importance(fm_dredge)

# Get top-most models, but fitted by REML (subset only models with SMW < 2)
(dd <- subset(fm_dredge, delta < 2))

# Generate confidence intervals
confset.95p <- get.models(fm_dredge, subset=delta < 2)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
confint(avgmod.95p)

####################################################################################
########## Prepare habitat suitability map (aka species distribution map) ##########
####################################################################################

# Load packages (if not installed, install from the repository)
if(!require(raster)) install.packages('raster',repos="http://cran.us.r-project.org");require(raster)
if(!require(rgdal)) install.packages('rgdal',repos="http://cran.us.r-project.org");require(rgdal)
if(!require(sp)) install.packages('sp',repos="http://cran.us.r-project.org");require(sp)
if(!require(rgeos)) install.packages('rgeos',repos="http://cran.us.r-project.org");require(rgeos)
if(!require(maps)) install.packages('maps',repos="http://cran.us.r-project.org");require(maps)
if(!require(GISTools)) install.packages('GISTools',repos="http://cran.us.r-project.org");require(GISTools)

# Import standardised raster files (these are predictor raster files)
ls() # I have packaged standardised raster files into an R data, so that they are all in one place.

# Alternatively, you can load individual file and standardise as follows:
# LC_CWED_4km <-raster("landcoverAll_CWED_4km.tif") # import raw raster
# LC_CWED_4kmstd <- (LC_CWED_4km-mean(det.data$LC_CWED_4km))/sd(det.data$LC_CWED_4km) # standardise for mapping

# Use coefficients from the top model to predict the distribution
logit.occ <- -1.33319 + (-0.17281*AGRI_GYRATE_1KMstd) + (0.63850*mean(scale(det.data$effort))) + (-0.13564*LC_CWED_4kmstd) + 
  (-0.25026*rivFM_2KMstd) + (-0.51087*roaFM_1KMstd) + (0.20839*sloposFM_4KMstd) + (0.21053*TC_CL2_GYRATE_8kmstd) + (0.13883*LC_GYRATE_1kmstd) +
  (-0.12507*BARE_GYRATE_2KMstd) + (0.06207*BLF_PD_32kmstd) + (-0.05952*BUILTUP_GYR_4kmstd)

# Backtransform (convert logit scale to probability scale)
occ<-exp(logit.occ)/(1+exp(logit.occ))

# Do the plots
plot(occ, col=rev(terrain.colors(100)))
plot(occ, col=rev(heat.colors(100)))

# Manual colour
pal <- colorRampPalette(c("blue", "cyan", "yellow", "orange", "red", "darkred"))
plot(occ, col=pal(100))

# Add boundary and camera trap locations (you can add capture locations but I don't recommend if the species is sensitive)
plot(Bhutan, add=T)
plot(sites, cex=0.4, col='grey20', add=T)

# Save the raster output (you can do post processing in ArcGIS/QGIS with this output - very helpful for management)
writeRaster(occ, "dhole_hab_new_jul19_SAC", format="GTiff")

###############################################
########## Model prediction accuracy ##########
###############################################

# Check the accuracy of your model
# AUC (area under the receiving operator characteristic curve) for model-averaged prediction
library(rms)
library(PresenceAbsence)
library(Hmisc)
library(HSAUR)
library(mgcv)
library(lattice)

# Load model-averaged raster output (map) and coordinates
mod <- occ
df <- det.data
coords <- df[, c(2, 3)]

# Add new columns to the data frame with prediction values for each site and assign ID
df$p <- extract(mod, coords)
df$ID <- seq.int(nrow(df))
data <- as.data.frame(cbind(df$ID, df$dhole, df$p))
names(data) <- c("id", "observed", "fitted")

# Check the accuracy and plot the results
presence.absence.hist(data, legend.cex=1, N.bars=10, opt.legend.cex=0.6, opt.thresholds=T, opt.methods=c(1))

accu1 <- presence.absence.accuracy(data, threshold=11, st.dev=F)
accu1[,-c(1,2)] <- signif(accu1[,-c(1,2)], digits=2)
accu1

accu2 <- presence.absence.accuracy(data, threshold=0.5, st.dev=F)
accu2[,-c(1,2)] <- signif(accu2[,-c(1,2)], digits=2)
accu2

# Do other plots for supplementary material
error.threshold.plot(data, color=T, main="Error Threshold Plot", opt.legend.cex=0.6, opt.methods=c(1, 4, 11))
auc.roc.plot(data, color=T, main="ROC Plot", opt.methods=c(1, 4, 11), opt.legend.cex=0.6)
optimal.thresholds(data, opt.methods=1:12)
calibration.plot(data, color=T)
#presence.absence.summary(data, opt.methods=c(1, 4, 11), N.bins=5, N.bars=10, legend.cex=0.6, opt.legend.cex=0.6)

## calculate optimal (kappa) threshold - this is helpful if you want to 
optimal.thresholds(data)

maxkappa <- as.numeric(optimal.thresholds(data, opt.methods="MaxKappa"))
( maxkappa <- maxkappa[2] )

( accu3 <- presence.absence.accuracy(data, threshold=maxkappa, st.dev=F) )

# Make a plot

library(visreg)
library(ggplot2)

visreg(fit=glm.mod1, xvar="elevation",  # use the best model here in place of glm.mod1
       scale="response", xlab="Elevation (standardised)", ylab="Habitat use prob") 

###################################################################################################################
####################################################### END #######################################################
###################################################################################################################

