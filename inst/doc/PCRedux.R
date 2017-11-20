## ---- echo=FALSE, fig.cap="\\label{samples_of_qPCRs}Sample data from A) the **C127EGHP** data set with 64 amplification curves (*chipPCR* package, [@rodiger_chippcr:_2015]) and B) from the **htPCR** data set with 8858 amplification curves (*qpcR* package, [@Ritz2008]).", fig.height=5.5, fig.width=11----
    
library(chipPCR)
suppressMessages(library(qpcR))


par(mfrow=c(1,2), las=0, bty="n", oma=c(1,1,1,1))
colors <- rainbow(ncol(C127EGHP)-2, alpha=0.3)
matplot(C127EGHP[, 2], C127EGHP[, c(-1, -2)], xlab="Cycle", ylab="RFU", 
        main="C127EGHP data set", type="l", lty=1, lwd=2, col=colors)
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

colors <- rainbow(ncol(htPCR)-1, alpha=0.2)
matplot(htPCR[, 1], htPCR[, c(-1)], xlab="Cycle", ylab="RFU", 
        main="htPCR data set", type="l", lty=1, lwd=2, col=colors)
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

## ---- echo=TRUE, fig.cap="\\label{autocorrelation}Autocorrelation analysis for amplification curves of the **testdat** data set (*qpcR* package, [@Ritz2008]). A) All amplification curves of the **testdat** data set. B) Positive curves and negative curves as determined by the `autocorrelation_test()` and a human. C) Performance analysis by the `performeR()` function (see Section \\ref{section_performeR} for details).", fig.height=5.5, fig.width=11----
# Test for autocorrelation in amplification curve data
# Load the libraries magrittr for pipes and qpcR for the data
library(magrittr)
library(qpcR)
library(PCRedux)
# Test for autocorrelation in the testdat data set
res_ac <- sapply(2L:ncol(testdat), function(i) {
    autocorrelation_test(testdat[, i])
}
)

# Plot curve data as overview
# Define the colors for the amplification curves
colors <- rainbow(ncol(testdat)-1, alpha=0.3)
# Names of samples
samples <- colnames(testdat)[-1]
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
matplot(testdat[, 1], testdat[, -1], xlab="Cycle", ylab="RFU", 
        main="testdat data set", type="l", lty=1, col=colors, lwd=2)
legend("topleft", samples, pch=19, col=colors, ncol=2, bty="n")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

# Curves rated by a human after analysis of the overview. 1 = positive, 
# 0 = negative
human_rating <- c(1,1,0,0,1,1,0,0,
                  1,1,0,0,1,1,0,0,
                  1,1,0,0,1,1,0,0)

# Convert the n.s. (not significant) in 0 and others to 1. 
# Combine the results of the aromatic autocorrelation_test as variable "ac", 
# the human rated values as variable "hr" in a new data frame (res_ac_hr).
res_ac_hr <- data.frame(ac=ifelse(res_ac=="n.s.", 0, 1), 
                        hr=human_rating) %>% as.matrix
                        res_performeR <- performeR(res_ac_hr[, "ac"], res_ac_hr[, "hr"])
                        
                        # Add ratings by human and autocorrelation_test to the plot
                        par(las=2)
                        plot(1:nrow(res_ac_hr), res_ac_hr[, "hr"], xlab="Sample", 
			     ylab="Decisions", xaxt="n", yaxt="n", pch=19)
                        axis(2, at=c(0,1), labels=c("negative", "positive"), las=2)
                        axis(1, at=1:nrow(res_ac_hr), labels=colnames(testdat)[-1], las=2)
                        points(1:nrow(res_ac_hr), res_ac_hr[, "ac"], pch=1, cex=2, 
			      col="red")
                        legend("topleft", c("Human", "autocorrelation_test"), pch=c(19,1), 
                               bty="n", col=c("black","red"))
                        mtext("B", cex = 2, side = 3, adj = 0, font = 2, las=0)
                        
                        barplot(as.matrix(res_performeR[, c(1:10,12)]), yaxt="n", ylab="", 
                                main="Performance of autocorrelation_test")
                        axis(2, at=c(0,1), labels=c("0", "1"), las=2)
                        mtext("C", cex = 2, side = 3, adj = 0, font = 2, las=0)

## ---- echo=FALSE, fig.cap="\\label{htPCR_nap_1}Sample data from the **htPCR** data set (*qpcR* package, [@Ritz2008]) with negative (black), ambiguous (red) or positive (green) amplification curves. The negative amplification curve has no sigmoid curvature, just a strong positive trend. The ambiguous amplification curve is not a perfect sigmoid curve, exhibits a positive slope in the background phase (cylce 1 $\rightarrow$ 5) and a low amplitude. In contrast, the positive amplification curve has a sigmoid curvature, starting with a background phase (cylce 1 $\rightarrow$ 5), an exponential phase (cylce 6 $\rightarrow$ 25) and a plateau phase (cylce 26 $\rightarrow$ 35).", fig.height=5.5, fig.width=11----
suppressMessages(library(qpcR))
matplot(htPCR[, 1], htPCR[, c(552, 512, 616)], xlab="Cycle", ylab="RFU", 
        main="htPCR data set", type="l", lty=1, lwd=2)
legend("topleft", c(paste("negative ", colnames(htPCR)[552]), 
                    paste("ambiguos ", colnames(htPCR)[512]), 
                    paste("positive ", colnames(htPCR)[616])
                   ), col=1:3, pch=19, bty="n")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
library(PCRedux)
suppressMessages(library(data.table))
library(magrittr)
filename <- system.file("decision_res_htPCR.csv", package="PCRedux")
decision_res_htPCR <- fread(filename) %>% as.data.frame()
head(decision_res_htPCR)
tail(decision_res_htPCR)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Use decision_modus to go rowise through all human ratings

dec <- lapply(1L:nrow(decision_res_htPCR), function(i) {
decision_modus(decision_res_htPCR[i, 2:4])
}) %>% unlist

names(dec) <- decision_res_htPCR[, 1]

# Show statistic of the decisions
summary(dec)

## ---- echo=FALSE, fig.cap="\\label{htPCR_nap}Summary of all sample data from the **htPCR** data set (*qpcR* package, [@Ritz2008]) with negative (black), ambiguous (red) or positive (green) amplification curves.", fig.height=5.5, fig.width=11----

# Plot the Frequencies of the decisions
barplot(table(dec), xlab="Decision", ylab="Frequency", 
        main="Amplification curves rated by humans\n htPCR data set", col=c(2,3,1))

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Decisions for sample P01.W06
decision_modus(decision_res_htPCR[which(decision_res_htPCR[["sample"]] == "P01.W06"), 
               2:4], max_freq=FALSE)

# Decisions for sample P01.W13
decision_modus(decision_res_htPCR[which(decision_res_htPCR[["sample"]] == "P01.W13"), 
               2:4], max_freq=FALSE)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Load the human rated data sets decision_res_htPCR.csv and the reevaluation 
# decision_res_htPCR_reevaluation.csv thereof.
# Load the decision_res_htPCR.csv data set
filename <- system.file("decision_res_htPCR.csv", package="PCRedux")
decision_res_htPCR <- as.data.frame(fread(filename))

# Load the decision_res_htPCR_reevaluation.csv data set
filename <- system.file("decision_res_htPCR_reevaluation.csv", package="PCRedux")
decision_res_htPCR_reevaluation <- as.data.frame(fread(filename))

# Show the first five elements of both data sets
head(decision_res_htPCR)
head(decision_res_htPCR_reevaluation)

# Renanme the columns to give them the same pattern
pattern <- c("sample", "test.result.1", "test.result.2", "test.result.3", "conformity")
colnames(decision_res_htPCR_reevaluation) <- pattern
colnames(decision_res_htPCR) <- pattern

# Merge the two data sets based on the sample name (amplification curve).
decision_res_htPCR_merged <- merge(decision_res_htPCR[, -5], 
                                   decision_res_htPCR_reevaluation[, -5], 
                                   by="sample")

# Show the first five elements of the merged data sets
head(decision_res_htPCR_merged)


## ---- echo=TRUE, fig.cap="\\label{decision_plot}decision plot.", fig.height=5.5, fig.width=11----

par(mfrow=c(1,3))


data <- decision_res_htPCR_merged
decision_plot <- function(data, name="", rating_columns=2:4){
    dec <- lapply(1L:nrow(data), function(i) {
        decision_modus(data[i, c(rating_columns)])
    }) %>% unlist
    barplot(table(dec), xlab="Decision", ylab="Frequency", 
            main=paste("Amplification curves rated by humans\n htPCR data set", 
                       name), col=c(2,3,1), ylim=c(0,6000))
}

decision_plot(decision_res_htPCR, name="", rating_columns=2:4)
decision_plot(decision_res_htPCR_reevaluation, name="", rating_columns=2:4)
decision_plot(decision_res_htPCR_merged, name="merged", rating_columns=2:7)

## ---- echo=TRUE, fig.cap="\\label{curve_fit_fail}Sample data from the **testdat** data set (*qpcR* package, [@Ritz2008]) with negative and positive amplification curves.", fig.height=7, fig.width=7----
# Load the drc package for the five-parameter log-logistic fit.
library(drc)
library(chipPCR)

# Load the testdat data set from the qpcR package without 
# loading the entire package.
data(testdat, package="qpcR")

# Arrange graphs in an orthogonal matrix and set the plot parmeters.
par(mfrow = c(2,2))

# Apply the the amptester function to all amplification curve data
# and write the results to the main of the plots.

curve_data_columns <- c(2,4,9,20)

for(i in 1L:length(curve_data_columns)) {
    res.ampt <- suppressMessages(amptester(testdat[, curve_data_columns[i]]))
    
    # Make a logical connection by two tests (shap.noisy, lrt.test and 
    # tht.dec) of amptester to decide if an amplification reaction is 
    # positive or negative.
    decision <- ifelse(!res.ampt@decisions[1] &&
                       res.ampt@decisions[2] && 
                       res.ampt@decisions[4], 
                       "positive", "negative")
    
    plot(testdat[, 1], testdat[, curve_data_columns[i]], 
         xlab = "Cycle", ylab = "RFU", pch = 19, main = "")
    mtext(LETTERS[i], cex = 2, side = 3, adj = 0, font = 2)
    legend("bottomright", paste0(colnames(testdat)[curve_data_columns[i]], 
                                 "\nDecision: ", decision), 
           bty="n", cex=1.25, col="red")
    # Use the drm function with a five-parameter log-logistic fit model.

        try(lines(predict(drm(testdat[, curve_data_columns[i]] ~ testdat[, 1], 
                          fct = LL.3())), col = 2), silent=TRUE)

        try(lines(predict(drm(testdat[, curve_data_columns[i]] ~ testdat[, 1], 
                          fct = LL.4())), col = 3), silent=TRUE)
 }

## ---- echo=TRUE, fig.cap="\\label{earlyreg_slopes}Amplification curves from the **htPCR** data set (*qpcR* package). The amplification curves show different slopes in the early phase (cycle 1 to 10) of the qPCR.", fig.height=5.5, fig.width=11----
par(bty="n", font=2, font.axis=2, las=1, cex.axis=1.3, cex.lab=1.3, lwd=1.3, 
    bg="white", oma=c(.5,.5,.5,.5))
data <- htPCR[, 1:501]
curve_colors <- c(rainbow(ncol(data)-1, alpha=.5))
range <- 1:10

par(mfrow=c(1,2), las=0, bty="n", oma=c(1,1,1,1))

matplot(data[, 1], data[, -1], col= curve_colors, pch=19, lty=1, type="l", 
        xlab="Cycle", ylab="RFU", main="htPCR data set")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

matplot(data[range, 1], data[range, -1], col= curve_colors, pch=19, lty=1, 
        type="l", xlab="Cycle", ylab="RFU", 
        main="htPCR data set\n Cycle 1 to 10")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

## ---- echo=TRUE, fig.cap="\\label{earlyreg_results}Clusters of samples according to their slope and intercept. The **htPCR** data set (*qpcR* package) was used. The amplification curves show different slopes in the early phase (cycle 1 to 10) of the qPCR. Both the slope and the intercept were used for a cluster analysis (k-means).", fig.height=7, fig.width=7----
# Normalize each amplification curve to their 0.99 percentile and use the 
# earlyreg function to determine the slope and intercept of the first 
# 6 cycles

res_slope_intercept <- do.call(rbind, lapply(2L:ncol(data), function(i){
                     earlyreg(x=data[, 1], y=data[, i], range=7, normalize=TRUE)
                    }))
# Label the sample with their original names
rownames(res_slope_intercept) <- colnames(htPCR)[2L:ncol(data)]
# Use k-means for cluster analysis
res_kmeans <- kmeans(res_slope_intercept[, "slope"], 3)$"cluster"

par(mfrow=c(2,3))
plot(res_slope_intercept, col=res_kmeans)
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

plot(density(res_slope_intercept[, "intercept"]), main="Intercept")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

plot(density(res_slope_intercept[, "slope"]), main="Slope")

for(i in unique(res_kmeans)) {
    abline(v=min(res_slope_intercept[res_kmeans==i, "slope"]), col=i)
}
mtext("C", cex = 2, side = 3, adj = 0, font = 2)

matplot(data[, 1], data[, which(res_kmeans==1)+1], col= curve_colors, pch=19, lty=1, 
 type="l", xlab="Cycle", ylab="RFU", 
 main="htPCR data set (1-500)\n Cluster 1", ylim=c(min(data[, -1]), 0.4))
mtext("D", cex = 2, side = 3, adj = 0, font = 2)

matplot(data[, 1], data[, which(res_kmeans==2)+1], col= curve_colors, pch=19, lty=1, 
 type="l", xlab="Cycle", ylab="RFU", 
 main="htPCR data set (1-500)\n Cluster 2", ylim=c(min(data[, -1]), 0.4))
mtext("E", cex = 2, side = 3, adj = 0, font = 2)

matplot(data[, 1], data[, which(res_kmeans==3)+1], col= curve_colors, pch=19, lty=1, 
 type="l", xlab="Cycle", ylab="RFU", 
 main="htPCR data set (1-500)\n Cluster 3", ylim=c(min(data[, -1]), 0.4))
mtext("F", cex = 2, side = 3, adj = 0, font = 2)

## ---- echo=TRUE, fig.cap="\\label{dil4reps94_head2tailratio}Calculation of the ratio between the head and the tail of a quantitative PCR amplification curve.", fig.height=5.5, fig.width=11----
suppressMessages(library(qpcR))


colors <- c("#FF00004D", "#80FF004D", "#00FFFF4D", "#8000FF4D")

data <- dil4reps94[, c(1, 2, 98, 194,290)]

res_head2tailratio <- lapply(2L:ncol(data), function(i) {
    head2tailratio(y=data[, i], normalize=TRUE, slope_normalizer=TRUE, 
                   verbose=TRUE)
})

data_normalized <- cbind(data[, 1], 
              sapply(2L:ncol(data), function(i){
                  data[, i] / quantile(data[, i], 0.999)
            })
             )

par(mfrow=c(1,2), las=0, bty="n", oma=c(1,1,1,1))

matplot(data_normalized[, 1], data_normalized[, -1], 
        xlab="Cycle", ylab="normalized RFU", main="dil4reps94 data set\nsubset", 
        type="l", lty=1, lwd=2, col=colors
       )
for(i in 1L:(ncol(data_normalized)-1)) {
    points(res_head2tailratio[[i]]$x_roi, res_head2tailratio[[i]]$y_roi, 
           col=colors[i], pch=19, cex=1.5)
    abline(res_head2tailratio[[i]]$fit, col=colors[i], lwd=2)
}

mtext("A", cex = 2, side = 3, adj = 0, font = 2)

colors <- c(rep("#FF00004D", 94),
            rep("#80FF004D", 94),
            rep("#00FFFF4D", 94),
            rep("#8000FF4D", 94)
            )

matplot(dil4reps94[, 1], dil4reps94[, -1], xlab="Cycle", ylab="RFU", 
        main="dil4reps94 data set", type="l", lty=1, lwd=2, col=colors)
mtext("B", cex = 2, side = 3, adj = 0, font = 2)


## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Calculate slope and intercept on noise (negative) amplification curve data
# for the last eight cycles.
library(qpcR)
library(magrittr)

res_hook <- sapply(2L:ncol(boggy), function(i) {
    hookreg(x=boggy[, 1], y=boggy[, i])}) %>% t %>% 
    data.frame(sample=colnames(boggy)[-1],.)


    data_colors <- rainbow(ncol(boggy[, -1]), alpha=0.5)
    cl <- kmeans(na.omit(res_hook[, 2:3]), 2)$cluster

## ---- eval=TRUE, echo=TRUE, results='asis'-------------------------------
library(xtable)
print(xtable(res_hook, caption = "Screening results for the analysis with the 
             hookreg algorithm. Samples withe a value of 1 in the hook column 
             had all a hook effect like curvature. The samples F4.1, F4.2, F5.1, 
             F5.2, F6.1 and F6.2 miss entries because the hoogreg algorithm 
             could not fit a linear model. This is an expected behavior, since these 
             amplification curve did not have a hook effect like curvature.", 
             label='tablehookreg'), include.rownames = FALSE)

## ---- echo=TRUE, fig.cap="\\label{plot_hookreg}Analysis of amplification curves for hook effect-like structures. ", fig.height=5, fig.width=7----
    par(mfrow=c(1,2))
    matplot(x=boggy[, 1], y=boggy[, -1], xlab="Cycle", ylab="RFU", 
            main="boggy Data Set", type="l", lty=1, lwd=2, col=data_colors)
    legend("topleft", as.character(res_hook$sample), pch=19, 
           col=data_colors, bty="n")
    mtext("A", cex = 2, side = 3, adj = 0, font = 2)
    
    plot(res_hook$intercept, res_hook$slope, pch=19, cex=2, col=data_colors,
         xlab="intercept", ylab="Slope", 
         main="Hook Effect-like Curvature\nboggy Data Set")
    mtext("B", cex = 2, side = 3, adj = 0, font = 2)
    points(res_hook$intercept, res_hook$slope, col=cl, pch=cl, cex=cl)
    legend("topright", c("Strong Hook effect", " Weak Hook effect"), 
           pch=c(1,2), col=c(1,2), bty="n")
    text(res_hook$intercept, res_hook$slope, res_hook$sample)

## ---- echo=TRUE, fig.cap="\\label{plot_mblrr}Robust local regression to analyze amplification curves. The amplification curves were arbitrarily selected from the htPCR data set. Not the differences in slopes and intercepts (red and green lines). The mblrr() function is presumably useful for data sets which are accompanied by noise and artifacts.", fig.height=7, fig.width=7----
suppressMessages(library(qpcR))
par(mfrow=c(3,2))
data <- htPCR[, c(1, 20, 500, 3000, 6000, 6500, 8000)]
for(i in 2L:ncol(data)){
        x <- data[, 1]
        y_tmp <- data[, i]/quantile(data[, i], 0.999)
        res_q25 <- y_tmp < quantile(y_tmp, 0.25)
        res_q75 <- y_tmp > quantile(y_tmp, 0.75)
        res_q25_lm <- try(suppressWarnings(lmrob(y_tmp[res_q25] ~ x[res_q25])), 
                          silent=TRUE)
        res_q75_lm <- try(suppressWarnings(lmrob(y_tmp[res_q75] ~ x[res_q75])), 
                          silent=TRUE)
        
        plot(x, y_tmp, xlab="Cylce", ylab="RFU (normalized)", 
             main=colnames(data)[i], 
             type="b", pch=19)
        abline(res_q25_lm, col="red")
        points(x[res_q25], y_tmp[res_q25], cex=2.5, col="red")
        abline(res_q75_lm, col="green")
        points(x[res_q75], y_tmp[res_q75], cex=2.5, col="green")
    }

## ---- eval=TRUE, echo=TRUE, results='asis'-------------------------------
# Load the library xtable for an appealing table output
library(xtable)

# Analyze the data via the mblrr() function

res_mblrr <- do.call(cbind, lapply(2L:ncol(data), function(i) {
                suppressMessages(mblrr(x=data[, 1], y=data[, i], 
                      normalize=TRUE)) %>% data.frame
}))
colnames(res_mblrr) <- colnames(data)[-1]

# Transform the data for a tabular output and assign the results to the object
# output_res_mblrr.

output_res_mblrr <- res_mblrr %>% t

# The output variable names of the mblrr() function are rather long. For better
# readability the variable names were changed to "nBG" (intercept of head region), 
# "mBG" (slope of head region), "rBG" (Pearson correlation of head region), 
# "nTP" (intercept of tail region), "mTP" (slope of tail region), "rBG" (Pearson 
# correlation of tail region)

colnames(output_res_mblrr) <- c("nBG", "mBG", "rBG",
                                "nTP", "mTP", "rTP")

print(xtable(output_res_mblrr, caption = "mblrr text intro. nBG, intercept of 
             head region; mBG, slope of head region; rBG, Pearson 
             correlation of head region; nTP, intercept of tail region; mTP, 
             slope of tail region; rBG, Pearson correlation of tail region", 
             label='tablemblrrintroduction'))


## ---- echo=FALSE, fig.cap="\\label{plot_raw_data_FFTrees}Fast and Frugal Trees.", fig.height=5.5, fig.width=11----
# Load the qpcR package to get some data
suppressMessages(library(qpcR))

# Use testdat data set, which contains positive and negative amplification 
# curves. Note, that noise was added to the raw data.

data <- data.frame(testdat[, 1], testdat[, -1] + rnorm(nrow(testdat[, -1])*ncol(testdat[, -1]), 0, 0.2))

# Visualize the data
par(mfrow=c(1,2), las=0, bty="n", oma=c(1,1,1,1))
colors <- rainbow(ncol(data)-2, alpha=0.3)
matplot(data[, 1], data[, -1], xlab="Cycle", ylab="RFU", 
        main="testdat data set", type="l", lty=1, lwd=2, col=colors)
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

# Give classes to the amplification curves (positive = TRUE, negative = FALSE)

human_rater_classification <- c(TRUE, TRUE,
                                FALSE, FALSE,
                                TRUE, TRUE,
                                FALSE, FALSE,
                                TRUE, TRUE,
                                FALSE, FALSE,
                                TRUE, TRUE,
                                FALSE, FALSE,
                                TRUE, TRUE,
                                FALSE, FALSE,
                                TRUE, TRUE,
                                FALSE, FALSE)

# Plot the data to see if the classification by the human is correct

colors <- as.factor(human_rater_classification)
matplot(data[, 1], data[, -1], xlab="Cycle", ylab="RFU", 
        main="testdat data set", type="l", lty=1, lwd=2, col=colors)

mtext("B", cex = 2, side = 3, adj = 0, font = 2)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Load the library xtable for an appealing table output
suppressMessages(library(FFTrees))

# Analyze the testdat data via the mblrr() function

res_mblrr <- do.call(cbind, lapply(2L:ncol(data), function(i) {
    suppressMessages(mblrr(x=data[, 1], y=data[, i], 
                           normalize=TRUE)) %>% data.frame
}))
colnames(res_mblrr) <- colnames(data)[-1]

# Transform the data for a tabular output and assign the results to the object
# output_res_mblrr.

output_res_mblrr <- res_mblrr %>% t

# The output variable names of the mblrr() function are rather long. For better
# readability the variable names were changed to "nBG" (intercept of head region), 
# "mBG" (slope of head region), "rBG" (Pearson correlation of head region), 
# "nTP" (intercept of tail region), "mTP" (slope of tail region), "rBG" (Pearson 
# correlation of tail region)

colnames(output_res_mblrr) <- c("nBG", "mBG", "rBG",
                                "nTP", "mTP", "rTP")

output_res_mblrr <- data.frame(class=human_rater_classification,
                          output_res_mblrr
                         )

output_res_mblrr.fft <- suppressMessages(FFTrees(formula = class ~., 
                                    data = output_res_mblrr[, c(1,2,3,5,6)]))

## ---- data=output_res_mblrr, eval=TRUE, echo=TRUE, results='asis'--------
library(xtable)
print(xtable(output_res_mblrr, caption = "mblrr text. nBG, intercept of 
             head region; mBG, slope of head region; rBG, Pearson 
             correlation of head region; nTP, intercept of tail region; mTP, 
             slope of tail region; rBG, Pearson correlation of tail region", 
             label='tablemblrr'), include.rownames = FALSE)

## ---- echo=FALSE, fig.cap="\\label{plot_FFTrees}Fast and Frugal Trees. nBG, intercept of head region; mBG, slope of head region; rBG, Pearson correlation of head region; nTP, intercept of tail region; mTP, slope of tail region; rBG, Pearson correlation of tail region.", fig.height=11, fig.width=11----
plot(output_res_mblrr.fft, decision.lables = c("Positive", "Negative"))

## ---- echo=TRUE, fig.cap="\\label{diffQ2_slope}Analysis of the amplification curve via diffQ2.", fig.height=5, fig.width=11----
# Load example data (sample F6.1) from the testdat data set
library(qpcR)
library(magrittr)
# Load MBmca package to calculate the minimum and the maximum of the second 
# derivative

library(MBmca)
data <- testdat[, c(1,22)]

# Calculate the minimum and the maximum of the second 
# derivative and assign it to the object res_diffQ2

res_diffQ2 <- suppressMessages(diffQ2(data, plot=FALSE, fct=min))

# Build a linear model from der second derivative of res_diffQ2
res_diffQ2_lm <- lm(res_diffQ2[["yTm1.2.D2"]] ~ res_diffQ2[["xTm1.2.D2"]])

par(mfrow=c(1,2))

data %>% plot(., xlab="Cycle", ylab="RFU", main="F6.1 (testdat)", type="l", 
                   lty=1, lwd=2, col="black")
abline(v=res_diffQ2[["xTm1.2.D2"]], col="grey", lwd=2)
mtext("A", cex = 2, side = 3, adj = 0, font = 2)


plot(res_diffQ2[["xTm1.2.D2"]], res_diffQ2[["yTm1.2.D2"]], pch=19, cex=1.5,
     xlab="Cycle", ylab="dd(RFU)/d(T)", 
     main="minimum and the maximum\n of second derivative")
abline(res_diffQ2_lm, col="blue", lwd=2)
legend("bottomright", paste("Slope: ", coefficients(res_diffQ2_lm)[2]), bty="n")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  # Copy and paste the code to an R console to evaluate it
#  
#  library(parallel)
#  library(doParallel)
#  
#  pcrfit_parallel <- function(data, n_cores = 1) {
#      # Determine the number of available cores and register them
#      if(n_cores == "all")
#          n_cores <- detectCores()
#  
#          registerDoParallel(n_cores)
#  
#          # Prepare the data for further processing
#          # Normalize RFU values to the alpha percentile (0.999)
#          cycles <- data.frame(cycles=data[, 1])
#          data_RFU <- data.frame(data[, -1])
#          data_RFU_colnames <- colnames(data_RFU)
#          data_RFU <- sapply(1L:ncol(data_RFU), function(i) {
#              data_RFU[, i] / quantile(data_RFU[, i], 0.999, na.rm=TRUE)
#          })
#          colnames(data_RFU) <- data_RFU_colnames
#  
#          # just to shut RCHeck for NSE we define ith_cycle
#          ith_cycle <- 1
#  
#          run_res <- foreach::foreach(ith_cycle = 1L:ncol(data_RFU),
#                             .packages=c("bcp", "changepoint", "chipPCR", "ecp", "MBmca",
#                                         "PCRedux", "pracma", "qpcR", "robustbase",
#                                         "zoo"),
#                             .combine = rbind) %dopar% {
#                                 suppressMessages(pcrfit_single(data_RFU[, ith_cycle]))
#                             }
#  
#  
#                             res <- cbind(runs = colnames(data_RFU), run_res)
#  
#                             rownames(res) <- NULL
#  
#                             res
#  
#  }
#  
#  # Calculate curve features of an amplification curve data. Note that not all
#  # available CPU cores are used. If need set "all" to use all available cores.
#  # In this example the testdat data set from the qpcR package is used.
#  # The samples F1.1 and F1.2 are positive amplification curves. The samples
#  # F1.3 and F1.4 are negative.
#  
#  library(qpcR)
#  res_pcrfit_parallel <- pcrfit_parallel(testdat[, 1:5])
#  res_pcrfit_parallel

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
# Calculate slope and intercept on noise (negative) amplification curve data.
# Load additional packages for data and pipes.
library(qpcR)
library(chipPCR)
library(fda.usc)
library(magrittr)

# Convert the qPCR data set to the fdata format
# Use unprocessed data from the testdat data set
res_fdata <- qPCR2fdata(testdat)

# Use preprocessed data (preprocess=TRUE) from the testdat data set
res_fdata_preprocessed <- qPCR2fdata(testdat, preprocess=TRUE)

# Extract column names and create rainbow color to label the data
res_fdata_colnames <- testdat[-1] %>% colnames()
data_colors <- rainbow(length(res_fdata_colnames), alpha=0.5)

## ---- echo=TRUE, fig.cap="\\label{qPCR2fdata}Clustering of amplification curves via Hausdorff distance. The amplification curves were preprocessed with the qPCR2fdata() function and subsequent processed by a cluster analysis and Hausdorff distance analysis.", fig.height=11, fig.width=11----
# Plot the converted qPCR data
par(mfrow=c(2,2))
res_fdata %>% plot(., xlab="Cycle", ylab="RFU", main="testdat", type="l", 
                   lty=1, lwd=2, col=data_colors)
legend("topleft", as.character(res_fdata_colnames), pch=19, 
       col=data_colors, bty="n", ncol=2)
mtext("A", cex = 2, side = 3, adj = 0, font = 2)
# Calculate the Hausdorff distance (fda.usc) package and plot the distances
# as clustered data.

res_fdata_hclust <- metric.hausdorff(res_fdata)
plot(hclust(as.dist(res_fdata_hclust)), main="Clusters of the amplification\n
curves as calculated by the Hausdorff distance")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)
rect(0.5,-3,12.25,0.5, border = "red")
text(7, 1, "negative", col="red")
rect(12.5,-3,24.5,0.5, border = "green")
text(14, 1, "positive", col="green")

# Repeat the plot and the cluster analysis for the preprocessed data
res_fdata_preprocessed %>% plot(., xlab="Cycle", ylab="RFU", main="testdat", type="l", 
                   lty=1, lwd=2, col=data_colors)
legend("topleft", as.character(res_fdata_colnames), pch=19, 
       col=data_colors, bty="n", ncol=2)
mtext("C", cex = 2, side = 3, adj = 0, font = 2)
# Calculate the Hausdorff distance (fda.usc) package and plot the distances
# as clustered data from preprocessed amplification curves.

res_fdata_hclust_preprocessed <- metric.hausdorff(res_fdata_preprocessed)
plot(hclust(as.dist(res_fdata_hclust_preprocessed)), 
     main="Clusters of the preprocessed amplification\n
curves as calculated by the Hausdorff distance")
rect(0.5,-3,12.25,0.5, border = "red")
text(7, 1, "negative", col="red")
rect(12.5,-3,24.5,0.5, border = "green")
text(14, 1, "positive", col="green")
mtext("D", cex = 2, side = 3, adj = 0, font = 2)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
# Calculate slope and intercept on positive amplification curve data from the
# VideoScan 32 cavity real-time PCR device.
# Load additional packages for data and pipes.
suppressMessages(library(data.table))
library(fda.usc)
library(magrittr)

# Load the qPCR data from the HCU32_aggR.csv data set

filename <- system.file("HCU32_aggR.csv", package="PCRedux")
data_32HCU <- fread(filename) %>% as.data.frame()

# Convert the qPCR data set to the fdata format
res_fdata <- qPCR2fdata(data_32HCU)

# Extract column names and create rainbow color to label the data
res_fdata_colnames <- data_32HCU[-1] %>% colnames()
data_colors <- rainbow(length(res_fdata_colnames), alpha=0.55)

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  # Load the qpcR package to calculate the Cq values by the second derivative
#  # maximum method.
#  
#  suppressMessages(library(qpcR))
#  
#  res_Cq <- sapply(2L:ncol(data_32HCU), function(i){
#      efficiency(pcrfit(data_32HCU, cyc=1, fluo=i, model=l6))
#  })
#  
#  data.frame(sample=colnames(data_32HCU)[-1], Cq=unlist(res_Cq["cpD2", ]), eff=unlist(res_Cq["eff", ]))
#  
#  #        Result
#  #
#  # sample    Cq      eff
#  # 1      A1 14.89 1.092963
#  # 2      B1 15.68 1.110480
#  # 3      C1 15.63 1.111474
#  # 4      D1 15.50 1.103089
#  # 5      E1 15.54 1.100122
#  # 6      F1 15.37 1.091367
#  # 7      G1 15.78 1.118713
#  # 8      H1 15.24 1.080062
#  # 9      A2 15.94 1.095004
#  # 10     B2 15.88 1.107878
#  # 11     C2 15.91 1.112694
#  # 12     D2 15.77 1.106286
#  # 13     E2 15.78 1.108201
#  # 14     F2 15.74 1.110697
#  # 15     G2 15.84 1.110749
#  # 16     H2 15.78 1.107229
#  # 17     A3 15.64 1.107543
#  # 18     B3 15.61 1.100984
#  # 19     C3 15.66 1.110703
#  # 20     D3 15.63 1.115996
#  # 21     E3 15.77 1.113885
#  # 22     F3 15.71 1.113985
#  # 23     G3 15.70 1.094108
#  # 24     H3 15.79 1.124223
#  # 25     A4 15.80 1.119774
#  # 26     B4 15.72 1.112124
#  # 27     C4 15.70 1.121453
#  # 28     D4 15.82 1.121809
#  # 29     E4 15.62 1.108028
#  # 30     F4 15.71 1.109634
#  # 31     G4 15.70 1.110373
#  # 32     H4 15.73 1.117827

## ---- echo=TRUE, fig.cap="\\label{HCU32}Clustering of amplification curves. The amplification curves from the 32HCU were processed with the ``qPCR2fdata()`` function and subsequent processed by a cluster analysis and Hausdorff distance analysis. A) Amplification curves were plotted from the raw data. B) Overall, the signal to noise ratios of the amplification curves were comparable between all cavities. C) The Cqs (Second Derivative Maximum) and the amplification efficiency (eff) were calculated with the ``efficiency(pcrfit())`` functions from the ``qpcR`` package. The median Cq is indicated as vertical line. Cqs larger or less than 0.1 of the Cq $\\tilde{x}$ are indicated with the labels of the corresponding sample. D) The clusters according to the Hausdorff distance show no specific pattern regarding the amplification curve signals. It appears that the samples D1, E1, F1, F3, G3 and H1 deviate most from the other amplification curves.", fig.height=8.5, fig.width=11----

library(fda.usc)
library(magrittr)

calculated_Cqs <- c(14.89, 15.68, 15.63, 15.5, 15.54, 15.37, 15.78, 15.24, 15.94, 
                    15.88, 15.91, 15.77, 15.78, 15.74, 15.84, 15.78, 15.64, 15.61, 
                    15.66, 15.63, 15.77, 15.71, 15.7, 15.79, 15.8, 15.72, 15.7, 15.82, 
                    15.62, 15.71, 15.7, 15.73)

calculated_effs <- c(1.09296326515231, 1.11047987547324, 1.11147389307153, 1.10308929700635, 
                     1.10012176315852, 1.09136717687619, 1.11871308210321, 1.08006168654712, 
                     1.09500422011318, 1.1078777171126, 1.11269436700649, 1.10628580163733, 
                     1.1082009954558, 1.11069683827291, 1.11074914659374, 1.10722949813473, 
                     1.10754282514113, 1.10098387264025, 1.1107026749644, 1.11599641663658, 
                     1.11388510347017, 1.11398547396991, 1.09410798249025, 1.12422338092929, 
                     1.11977386646464, 1.11212436173214, 1.12145338871426, 1.12180879952503, 
                     1.1080276005651, 1.10963449004393, 1.11037302758388, 1.11782689816295)

# Plot the converted qPCR data
layout(matrix(c(1,2,3,4,4,4), 2, 3, byrow = TRUE))
res_fdata %>% plot(., xlab="Cycle", ylab="RFU", main="HCU32_aggR", type="l", 
                   lty=1, lwd=2, col=data_colors)
legend("topleft", as.character(res_fdata_colnames), pch=19, 
       col=data_colors, bty="n", ncol=4)
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

# Plot the background and plateau phase.

boxplot(data_32HCU[,-1] - apply(data_32HCU[, -1], 2, min), 
        col=data_colors, las=2, main="Signal to noise ratio", 
        xlab="Sample", ylab="RFU")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

# Plot the Cqs and the amplification efficiencies.
# Determine the median of the Cq values and label all Cqs, which a less 0.1 Cqs 
# of the median or more then 0.1 Cqs of the median Cq.

plot(calculated_Cqs, calculated_effs, xlab="Cq (SDM)", 
     ylab="eff", main="Cq vs. Amplification Efficiency", 
     type="p", pch=19, lty=1, lwd=2, col=data_colors)

median_Cq <- median(calculated_Cqs)
abline(v=median_Cq)

text(median_Cq + 0.01, 1.085, expression(paste(tilde(x))))
labeled <- c(which(calculated_Cqs < median_Cq - 0.1),  
             which(calculated_Cqs > median_Cq + 0.1))

text(calculated_Cqs[labeled], calculated_effs[labeled], as.character(res_fdata_colnames)[labeled])
mtext("C", cex = 2, side = 3, adj = 0, font = 2)

# Calculate the Hausdorff distance using the fda.usc package and cluster the 
# the distances.

res_fdata_hclust <- metric.hausdorff(res_fdata)
cluster <- hclust(as.dist(res_fdata_hclust))

# plot the distances as clustered data and label the leafs with the Cq values and
# colored dots.

plot(cluster, main="Clusters of the amplification\n
curves as calculated by the Hausdorff distance")
mtext("D", cex = 2, side = 3, adj = 0, font = 2)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Calculate curve features of an amplification curve data set.
# Use the C126EG685 data set from the chipPCR package and analyze the samples
# A01, A02, A04 and B05.

library(chipPCR)

res_1 <- cbind(runs="A01", pcrfit_single(C126EG685[, 2]))
res_2 <- cbind(runs="A02", pcrfit_single(C126EG685[, 3]))
res_3 <- cbind(runs="A04", pcrfit_single(C126EG685[, 5]))
res_4 <- cbind(runs="B04", pcrfit_single(C126EG685[, 17]))


res <- rbind(A01=res_1, A02=res_2, A04=res_3, B04=res_4)

## ---- echo=TRUE, fig.cap="\\label{visdat_pcrfit_plot}Application of visdat_pcrfit() for the visualization of the data structure after an analysis by pcrfit_single().", fig.height=11, fig.width=11----
# Show all results in a plot. Not that the interactive parameter is set to 
# FALSE.

visdat_pcrfit(res, type="all", interactive=FALSE)

