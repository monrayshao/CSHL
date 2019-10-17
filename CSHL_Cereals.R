# Author: Mon-Ray Shao
# Institute: Donald Danforth Plant Science Center
# Date: October 17, 2019
# Course: CSHL Cereal Genomics (Topp session)

# Description:
# Example R code to perform analyis of root traits from GIA Roots output



# NAM GIA ROOTS -----------------------------------------------------------
# Notes:
# A LemnaTec trial was used to grow seedlings from the Maize NAM founders
# Seedlings were watered under 25%, 50%, 75%, or 100% soil water capacity
# After 24 days, the seedlings were removed from their pots and roots washed
# Root systems were then imaged using a digital camera and processed with GIA roots
# The resulting data was exported and used for this demonstration



# Import Data and Inspect -------------------------------------------------

gia <- read.csv("https://raw.githubusercontent.com/monrayshao/CSHL/master/GIA_CSHL_Cereal.csv", stringsAsFactors = F)
head(gia) # Preview data frame
dim(gia) # Print dimensions of data frame
colnames(gia) # Print column names of data frame

gia$Network.Area # Extract a single column
mean(gia$Network.Area) # Taking the mean
mean(subset(gia, treatment == 50)$Network.Area) # One way to subset data
mean(gia[which(gia$treatment == 50), ]$Network.Area) # Another way to subset data
aggregate(Network.Area ~ genotype + treatment, gia, FUN = mean) # Calculate means according to multiple variables

class(gia$treatment) # R assumes the Treatment column is numeric
gia$treatment <- as.factor(gia$treatment) # To avoid confusion in regression vs classification, convert to factor

plot(gia[, 7:12]) # R easily makes scatterplots
plot(subset(gia, treatment == 100)[, 7:12]) # Combine with subsetting data

# Convert to long data format using tidyr, for ggplot2 plotting
install.packages("tidyr")
library(tidyr)
dataCols <- 7:44 # Save trait columns #'s for convenience
gatherData <- gather(gia, key = "Traits", value = "Value", dataCols) # Gather turns wide to long
head(gatherData); tail(gatherData)
class(gatherData$Value) # Check that all phenotypes converted were indeed numeric, avoid data class mixing

# Make boxplot of every trait with faceted windows
install.packages("ggplot2")
library(ggplot2)
ggplot(gatherData, aes(x = treatment, y = Value)) + 
  geom_boxplot(aes(color = treatment)) + 
  facet_wrap(~Traits, scale = "free_y")

# Or jitter plot
ggplot(gatherData, aes(x = treatment, y = Value)) + 
  geom_jitter(width = 0.1, height = 0, alpha = 0.25, aes(color = treatment)) + 
  facet_wrap(~Traits, scale = "free_y")



# Data Quality ------------------------------------------------------------

# Histograms
# Use a loop to make a density plot of every trait
par(mfrow = c(3,3)) # Make 3x3 plots per page
for (i in dataCols){
  gia25 <- na.omit(subset(gia, treatment == 25)[,i]) # Subset the treatments for each trait
  gia50 <- na.omit(subset(gia, treatment == 50)[,i])
  gia75 <- na.omit(subset(gia, treatment == 75)[,i])
  gia100 <- na.omit(subset(gia, treatment == 100)[,i])
  ymax <- max(density(gia25)$y, density(gia50)$y, density(gia75)$y, density(gia100)$y) # Establish y-axis for each trait
  
  plot(density(na.omit(gia[,i])), main = colnames(gia)[i], col = NA, ylim=c(0, ymax)) # Plot empty frame of right size
  polygon(density(gia25), col = rgb(1, 0, 0, 0.25)) # Now plot each density curve
  polygon(density(gia50), col = rgb(1, 1, 0, 0.25)) # With rgb() you can set transparency
  polygon(density(gia75), col = rgb(0, 1, 0, 0.25))
  polygon(density(gia100), col = rgb(0, 0, 1, 0.25))
}
par(mfrow = c(1,1)) # Go back to 1 plot per page

# QQ Plots
par(mfrow = c(3,3))
for (i in dataCols){
  gia25 <- na.omit(subset(gia, treatment == 25)[,i]) # Subset the treatments for each trait
  gia50 <- na.omit(subset(gia, treatment == 50)[,i])
  gia75 <- na.omit(subset(gia, treatment == 75)[,i])
  gia100 <- na.omit(subset(gia, treatment == 100)[,i])
  ymin <- min(qqnorm(gia25, plot.it = FALSE)$y, 
              qqnorm(gia50, plot.it = FALSE)$y,
              qqnorm(gia50, plot.it = FALSE)$y,
              qqnorm(gia100, plot.it = FALSE)$y)
  ymax <- max(qqnorm(gia25, plot.it = FALSE)$y, 
              qqnorm(gia50, plot.it = FALSE)$y,
              qqnorm(gia50, plot.it = FALSE)$y, 
              qqnorm(gia100, plot.it = FALSE)$y) # Establish y-axis for each trait
  
  qqnorm(na.omit(gia[,i]), main = colnames(gia)[i], col = NA, ylim = c(ymin, ymax)) # Plot empty frame of right size
  points(qqnorm(gia25, plot.it = FALSE), col = "red") # Now make each QQ plot
  points(qqnorm(gia50, plot.it = FALSE), col = "orange")
  points(qqnorm(gia75, plot.it = FALSE), col = "green")
  points(qqnorm(gia100, plot.it = FALSE), col = "blue")
  qqline(gia25, col = "red") 
  qqline(gia50, col = "orange")
  qqline(gia75, col = "green")
  qqline(gia100, col = "blue")
}
par(mfrow = c(1,1))

# Outliers, univariate example
subgia <- subset(gia, treatment == 25) # Just subset this as example
iqr <- IQR(subgia[, "Network.Area"]) # Just take Network.Area as example
perc25 <- quantile(subgia[,"Network.Area"])['25%'] # Calculate 1st quantile
perc75 <- quantile(subgia[,"Network.Area"])['75%'] # Calculate 3rd quantile
# Find outliers using the 1.5x IQR method
outlier <- which(subgia[,"Network.Area"] > perc75 + 1.5*iqr | subgia[,"Network.Area"] < perc25 - 1.5*iqr)
print(subgia[outlier,])

# Multivariate approach, one example
# Subset each treatment without missing data (function won't accept)
gia25 <- na.omit(subset(gia, treatment == 25))
gia50 <- na.omit(subset(gia, treatment == 50))
gia75 <- na.omit(subset(gia, treatment == 75))
gia100 <- na.omit(subset(gia, treatment == 100))
# Calculate Mahalanobis distance for each treatment
dist25 <- mahalanobis(gia25[ ,dataCols], center = colMeans(gia25[ ,dataCols]), cov = cov(gia25[ ,dataCols]), tol = 1e-20)
dist50 <- mahalanobis(gia50[ ,dataCols], center = colMeans(gia50[ ,dataCols]), cov = cov(gia50[ ,dataCols]), tol = 1e-20)
dist75 <- mahalanobis(gia75[ ,dataCols], center = colMeans(gia75[ ,dataCols]), cov = cov(gia75[ ,dataCols]), tol = 1e-20)
dist100 <- mahalanobis(gia100[ ,dataCols], center = colMeans(gia100[ ,dataCols]), cov = cov(gia100[ ,dataCols]), tol = 1e-20)
# Find which data points are beyond a given threshold
outliers25 <- which(dist25 > qchisq(0.999, df = length(dataCols)))
outliers50 <- which(dist50 > qchisq(0.999, df = length(dataCols)))
outliers75 <- which(dist75 > qchisq(0.999, df = length(dataCols)))
outliers100 <- which(dist100 > qchisq(0.999, df = length(dataCols)))
# Show these outliers, are they of the same genotype or not?
print(gia25[outliers25,]$genotype)
print(gia50[outliers50,]$genotype)
print(gia75[outliers75,]$genotype)
print(gia100[outliers100,]$genotype)
# Make list of outlier barcodes for removal
bad <- c(gia25[outliers25, ]$plantbarcode, 
         gia50[outliers50, ]$plantbarcode, 
         gia75[outliers75, ]$plantbarcode, 
         gia100[outliers100, ]$plantbarcode)
# Remove outliers
giaFiltered <- gia[-which(gia$plantbarcode %in% bad), ]
dim(giaFiltered)



# Trait Significance ------------------------------------------------------

# Correlations
install.packages("corrplot")
library(corrplot)
gia2 <- na.omit(giaFiltered) # Use no-missing data set for simplicity
corrmatrix <- cor(gia2[, dataCols], method = "spearman") # Calculate Spearman correlations
corrplot(corrmatrix, tl.cex = 0.5, tl.col = "black", cl.cex = 0.5) # Plot correlation matrix
corrplot(corrmatrix, tl.cex = 0.5, tl.col = "black", cl.cex = 0.5, order = "hclust", hclust.method = "complete")

# Heritability
install.packages("plyr")
install.packages("lme4")
library(plyr)
library(lme4)
# Make empty data frame to hold results
VarComp <- data.frame(Trait = colnames(giaFiltered)[dataCols],
                      G = NA, E = NA, GxE = NA, error = NA)
for (i in dataCols){
  form <- paste0( colnames(giaFiltered)[i], "~ 1 + (1|genotype) + (1|treatment) + (1|genotype:treatment)" )
  mod <- lmer( as.formula(form) , data = giaFiltered) # Run random effect model to extract variances
  # Error catcher below, we will remove traits whose models did not converge
  tryCatch( {lmer( as.formula(form) , data = giaFiltered)}, warning = function(w) print(colnames(giaFiltered)[i]))
  result <- as.data.frame(VarCorr(mod)) # Extract variance components
  varTotal <- sum(result$vcov)
  G <- result[ which(result$grp == "genotype"), ]$vcov
  E <- result[ which(result$grp == "treatment"), ]$vcov
  GxE <- result[ which(result$grp == "genotype:treatment"), ]$vcov
  error <- result[ which(result$grp == "Residual"), ]$vcov
  VarComp[i-6,]$G <- G / varTotal # Parse variance components into proper place in data frame
  VarComp[i-6,]$E <- E / varTotal
  VarComp[i-6,]$GxE <- GxE / varTotal
  VarComp[i-6,]$error <- error / varTotal
}  
questionable <- c("Maximum.Number.of.Roots", "Minor.Ellipse.Axis", "Network.Perimeter", "Specific.Root.Length_Excised")
VarComp <- VarComp[-which(VarComp$Trait %in% questionable), ] # Remove these traits whose models did not run properly
environments <- 4
replicates <- plyr::count(giaFiltered$genotype)
avgReplicates <- mean(replicates$freq)
# Calculate broad-sense heritability
VarComp$heritability <- VarComp$G / (VarComp$G + VarComp$GxE/environments + VarComp$error/avgReplicates)

# Re-ordering factor levels; can be a pain!
VarComp <- VarComp[rev(order(VarComp$heritability)),]
VarComp$Trait <- as.factor(VarComp$Trait)
ord <- match( VarComp[rev(order(VarComp$heritability)),]$Trait, levels(VarComp$Trait))
VarComp$Trait <- factor(VarComp$Trait, levels(VarComp$Trait)[ord])

install.packages("viridis")
library(viridis)
# Plot heritability estimates
ggplot(VarComp, aes(x = Trait, y = heritability)) +
  geom_bar(stat = "identity", color = "black", aes(fill = heritability)) + ylim(0,1) +
  labs(x = "", y = "Broad-Sense Heritability\n") + guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(option = "magma")

# Reverse levels to make nice horizontal plot
VarComp$Trait <- factor(VarComp$Trait, rev(levels(VarComp$Trait)))
ggplot(VarComp, aes(x = Trait, y = heritability)) +
  geom_bar(stat = "identity", color = "black", aes(fill = heritability)) + ylim(0,1) +
  labs(x = "", y = "\nBroad-Sense Heritability") + guides(fill = FALSE) +
  scale_fill_viridis(option = "magma") + coord_flip()

# Analysis of Variance for factor significance
# Make empty data frame
anovaRes <- data.frame(Trait = colnames(giaFiltered)[dataCols],
                      G = NA, E = NA, GxE = NA)
for (i in dataCols){
  form <- paste0( colnames(giaFiltered)[i], "~ genotype + treatment + genotype:treatment" ) 
  mod <- lm( as.formula(form) , data = giaFiltered) # Run linear model
  result <- as.data.frame(anova(mod))
  G <- result[ which(rownames(result) == "genotype"), ]$'Pr(>F)' # Extract ANOVA results for each factor
  E <- result[ which(rownames(result) == "treatment"), ]$'Pr(>F)'
  GxE <- result[ which(rownames(result) == "genotype:treatment"), ]$'Pr(>F)'
  anovaRes[i-6,]$G <- G # Parse results into correct place in data frame
  anovaRes[i-6,]$E <- E
  anovaRes[i-6,]$GxE <- GxE
}  
anovaRes$G_adj <- p.adjust(anovaRes$G, method = "fdr") # Adjust p-values!
anovaRes$E_adj <- p.adjust(anovaRes$E, method = "fdr")
anovaRes$GxE_adj <- p.adjust(anovaRes$GxE, method = "fdr")
gatherAnova <- gather(anovaRes, key = "Factor", value = "p.value", 5:7) # Convert to long table for ggplot
gatherAnova$log_p.value <- -log10(gatherAnova$p.value) # Plot p-values using -log10
# Plot heatmap of ANOVA results
ggplot(gatherAnova, aes(x = Trait, y = Factor)) +
  geom_tile(aes(fill = log_p.value)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top") +
  scale_fill_viridis(name = "-log10(p-value)")
  


# Principal Component Analysis --------------------------------------------

# Using prcomp()
gia.scaled <- scale(gia2[, dataCols], center = T, scale = T) # Since we have different types of features, scale the data
gia.pca <- prcomp(gia.scaled) # Run PCA function
screeplot(gia.pca, type = "lines") # Shows how much variance each PC explains
summary(gia.pca)$importance

# Parse the results for plotting
gia.pci <- data.frame(gia.pca$x)
gia.pci$genotype <- gia2$genotype
gia.pci$treatment <- gia2$treatment
gia.pci$treatment <- as.factor(gia.pci$treatment)

# Plot subset of PCA results
ggplot(subset(gia.pci, treatment == 25 | treatment == 100), 
       aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = treatment, color = treatment, fill = treatment)) +
  theme_classic()

# Plot full PCA results
ggplot(gia.pci, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = treatment, color = treatment, fill = treatment)) +
  stat_ellipse(type = "t", geom = "polygon", alpha = 0.25, aes(fill = treatment), color = NA) +
  labs(x = paste0("\nPC1 (", 100*summary(gia.pca)$importance[2,1], "% Variance)"),
       y = paste0("PC2 (", 100*summary(gia.pca)$importance[2,2], "% Variance)\n")) +
  scale_shape_manual(values = c(21 ,22, 23, 24)) +
  theme_classic()
# PC1 seems to separate treatments
# PC2 separates samples within each treatment -- possibly genotype?

# What traits contribute to each axis?
sort(gia.pca$rotation[,1]) # size-related traits
sort(gia.pca$rotation[,2]) # ratios, excised roots
biplot(gia.pca, cex = 0.75)


# Classification ----------------------------------------------------------

# Choose 80% of rows for training data
set.seed(100) # Set a seed so you get the same results each time you run on your computer
train <- sample(nrow(gia2), size = 0.8 * nrow(gia2), replace = FALSE) # Random observations for sample training set

# Random forest
install.packages("randomForest")
install.packages("rminer")
library(randomForest)
library(rminer)
rfData <- gia2[, c(5, dataCols)] # Only want treatment ID's and trait data
rfData$treatment <- as.factor(rfData$treatment) # Make sure treatment is data type "factor", not numeric!
trainData <- rfData[train,] # Separate training data
testData <- rfData[-train,] # Separate test data
set.seed(999) # Set seed again
rfModel <- randomForest(treatment ~ ., data = trainData, ntree = 5000, proximity = TRUE) # Run random forest
rfPredict <- predict(rfModel, testData) # Predict test data using random forest model
mmetric(testData$treatment, rfPredict, "ACC") # Check accuracy on test data
varImpPlot(rfModel, n.var = length(dataCols), sort = T) # Check variable importance

# Plot "visualization" of random forest
gia.mds <- data.frame(cmdscale(d = (1 - rfModel$proximity), k = 2)) # Calculate distances from proximity matrix
gia.mds$treatment <- as.factor(gia2[train, ]$treatment) # Parse treatment info
ggplot(gia.mds, aes(x = X1, y = X2)) +
  geom_point(aes(color = treatment, shape = treatment)) +
  stat_ellipse(type = "t", geom = "polygon", alpha = 0.25, aes(fill = treatment), color = NA) +
  theme_classic()



# Clustering --------------------------------------------------------------

# Hierarchical clustering example for Network.Area

area <- aggregate(Network.Area ~ genotype + treatment, gia2, FUN = mean) # Need to aggregate results for this
areaTable <- spread(area, key = "treatment", value = "Network.Area") # Convert to wide data for plotting
areaTable <- na.omit(areaTable) # Remove missing data for now
rownames(areaTable) <- areaTable$genotype # Add labels to rownames
areaTable <- as.matrix(areaTable[, -1]) # Don't want to plot this
d <- dist(areaTable, method = "euclidean") # Calculate distances
hc <- hclust(d, method = "complete") # Run agglomeration 
plot(hc) # Plot clustering results

install.packages("gplots")
library(gplots)
# Plot clustering heatmap (separates the treatments so we can see)
heatmap.2(t(areaTable), Rowv = FALSE, dendrogram = "column",
          margins = c(5, 7), trace = "none", col = "cm.colors", 
          main = "Network Area")




