## ---- echo=FALSE, fig.cap="\\label{figure_sigmoid_cuve} Amplification curve data from an iQ5 (Bio-Rad) thermo-cycler and a high throughput experiment in the Biomark HD (Fluidigm). A) The `C127EGHP` dataset (\\texttt{chipPCR} package, [@roediger2015chippcr]) with 64 amplification curves was produced in conventional thermo-cycler with a 8 x 12 PCR grid. B) The `htPCR` dataset (\\texttt{qpcR} package, [@Ritz2008]), which contains 8858 amplification curves, was produced in a 95 x 96 PCR grid. Only 200 amplification curves are shown. In contrast to `A)` have all amplification curves in `B)` an off-set (intercept) of circa 0.25 RFU. C) Model function of a one-parameter sigmoid function. D) Model function of a sigmoid function with an intercept $n$ = 0.2 RFU. E) Model function of a sigmoid function with an intercept ($n$ ~ 0.25 RFU) and a square portion $m * x^{2}$. F) Model function of a sigmoid function with an intercept ($n$) and a square portion of $m * x^{2}$ and additional noise $\\epsilon$ (normally distributed).", fig.height=8----
options(warn = -1)
library(chipPCR)
suppressMessages(library(qpcR))

x_val <- seq(-10, 10, 0.05)
y_val <- 1 / (1 + exp(-x_val))
y_val_slope <- 1 / (1 + exp(-x_val)) + 0.2
y_val_slope_quadratic <- 1 / (1 + exp(-x_val)) + 0.0005 * x_val ^ 2 + 0.2
y_val_slope_quadratic_noise <- 1 / (1 + exp(-x_val)) + 0.0005 * x_val ^ 2 + 0.2 + rnorm(length(x_val), mean = 0.01, sd = 0.05)

y_lim <- c(-0.05, max(c(
  y_val, y_val_slope, y_val_slope_quadratic,
  y_val_slope_quadratic_noise
)) * 1.2)

par(mfrow = c(3, 2), las = 0, bty = "o", oma = c(0, 0, 0, 0))

colors <- rainbow(ncol(C127EGHP) - 2, alpha = 0.5)
matplot(
  C127EGHP[, 2], C127EGHP[, c(-1, -2)], xlab = "Cycle", ylab = "RFU",
  main = "C127EGHP dataset", type = "l", lty = 1, lwd = 2, col = colors
)
abline(h = 0, col = "grey")
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)

colors <- rainbow(200, alpha = 0.5)
matplot(
  htPCR[, 1], htPCR[, c(2L:201)], xlab = "Cycle", ylab = "RFU", ylim = c(0, 1.75),
  main = "htPCR dataset", type = "l", lty = 1, lwd = 2, col = colors
)
abline(h = 0, col = "grey")
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)

plot(x_val, y_val, type = "l", xlab = "x", ylab = "f(x)", ylim = y_lim)
abline(h = 0, col = "grey")
legend("topleft", expression(y == frac(1, (1 + e ^ {
  -x
}))), bty = "n", cex = 0.9)
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2)

plot(x_val, y_val_slope, type = "l", xlab = "x", ylab = "f(x)", ylim = y_lim)
abline(h = 0, col = "grey")
legend("topleft", expression(y == frac(1, (1 + e ^ {
  -x
})) + n), bty = "n", cex = 0.9)
mtext("D", cex = 1.2, side = 3, adj = 0, font = 2)

plot(
  x_val, y_val_slope_quadratic, type = "l", xlab = "x", ylab = "f(x)",
  ylim = y_lim
)
abline(h = 0, col = "grey")
legend("topleft", expression(y == frac(1, (1 + e ^ {
  -x
})) + m * x ^ 2 + n), bty = "n", cex = 0.9)
mtext("E", cex = 1.2, side = 3, adj = 0, font = 2)

plot(
  x_val, y_val_slope_quadratic_noise, type = "l", xlab = "x", ylab = "f(x)",
  ylim = y_lim
)
abline(h = 0, col = "grey")
legend("topleft", expression(y == frac(1, (1 + e ^ {
  -x
})) + m * x ^ 2 + n + epsilon, epsilon %~% N(0, sigma)), bty = "n", cex = 0.9)
mtext("F", cex = 1.2, side = 3, adj = 0, font = 2)

## ---- echo=FALSE, fig.cap="Regions of interest in amplification curves. A) In general, the fluorescence emitted (RFU, relative fluoresce units) by the reporter dye (e.g, SYBR Green, EvaGreen) is plotted against cycle number. The amplification curve data was taken from the `testdat` dataset (\\texttt{qpcR} package). More generally, amplification curves can be divided into three regions of interest. These are the ground phase, exponential phase and plateau phase. `top`, takeoff point. `tdp`, takdown point. `sd_bg` is the standard deviation within the ground phase. The exponential region (red dots) can be used to determine the Cq values and estimates of amplification efficiency. The straight red line is the regression line of a linear model. In principle, after further processing steps (e.g., logarithmetic), slopes in this range can be determined. B) PCRs without amplification reaction are usually characterized by a flat (non-sigmoid) signal. C) The exponential phase of PCR reactions may vary considerably. Ideally, the slopes are the same for all reactions. This would be synonymous with the same amplification efficiency in all reactions. However, in practice, amplification curves with different increases are usually found. In particular, amplification curves that become detectable in later cycles often have lower increases.\\label{amplification_curve_ROI}", fig.height=8, fig.width=8----
options(warn = -1)
library(qpcR)
library(PCRedux)

colors <- rainbow(10, alpha = 0.15)

x_range <- 1L:35
d <- testdat[x_range, ]
amp_data <- data.frame(
  d[, 1],
  pos = d[, 3] + 0.9,
  posReverse = (max(d[, 3]) - rev(d[, 3])) + 0.9,
  neg = d[, 4] + 0.9 + 0.0005 * d[, 1] ^ 2
)

# Calculation for the normal data
res_amp_data <- pcrfit(amp_data, 1, 2, l5)
res_takeoff <- takeoff(res_amp_data)

# Calculation of sd_bg

res_sd_bg <- sd(amp_data[1:res_takeoff[[1]], 2])

# Calculation for the reversed data
res_amp_data_reverse <- pcrfit(amp_data, 1, 3, l5)
res_takeoff_reverse <- takeoff(res_amp_data_reverse)
res_takeoff_reverse[[1]] <- nrow(d) - res_takeoff_reverse[[1]]
res_takeoff_reverse[[2]] <- amp_data[res_takeoff_reverse[[1]], 2] - res_takeoff_reverse[[2]] + min(amp_data[, 3])

exponentialRange <- c((res_takeoff[[1]] + 1):(res_takeoff_reverse[[1]] - 1))

backgroundplateu <- function(x) {
  bg <- mean(head(x, res_takeoff[[1]])) + 3 * sd(head(x, res_takeoff[[1]]))
  plat <- mean(tail(x, 10)) - 3 * sd(tail(x, 10))
  list(bg = bg, plateau = plat)
}

res_lm <- lm(amp_data[exponentialRange, 2] ~ amp_data[exponentialRange, 1])

y_lim <- max(amp_data[, 2:4]) * 1.15

res_bgpl <- unlist(backgroundplateu(amp_data[, 2]))

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE), respect = TRUE)

plot(amp_data[, 1], amp_data[, 2], ylim = c(-0.1, y_lim), xlab = "Cycle", ylab = "RFU", type = "b", lwd = 2, pch = 19)

text(c(2,30), c(10.5,10.5), c("Head", "Tail"), cex = 1.2, col = "red")

rect(0, 0, res_takeoff[[1]] + 1, res_takeoff[[2]] * 1.25, col = colors[1], border = NA)
text(5, res_bgpl[1] * 1.45, "Ground phase")

rect(res_takeoff_reverse[[1]] - 1, res_takeoff_reverse[[2]] * 0.95, nrow(amp_data), y_lim, col = colors[5], border = NA)
text(32, res_bgpl[2] * 1.1, "Plateau phase")

text(res_takeoff_reverse[[1]], mean(amp_data[, 2]), "Exponential\nregion")


points(
  c(res_takeoff[[1]], res_takeoff_reverse[[1]]),
  c(res_takeoff[[2]], res_takeoff_reverse[[2]]), pch = 12, cex = 2.5
)

text(
    c(res_takeoff[[1]], res_takeoff_reverse[[1]]),
    c(res_takeoff[[2]], res_takeoff_reverse[[2]]) + c(1.05, -1.05), c("top", "tdp")
)


arrows(20, 0, 20, res_bgpl[1], code = 3, length = 0.1)
text(30, res_bgpl[1] / 2, "Background")

arrows(5, res_bgpl[2], 5, max(amp_data[, 2]), code = 3, length = 0.1)
text(15, res_bgpl[2] * 0.95, "Plateau")


abline(res_lm, col = "red")
points(amp_data[exponentialRange, 1], amp_data[exponentialRange, 2], pch = 19, col = "red")

abline(h = res_bgpl, col = c("green", "blue"))
abline(h = 0, col = "grey")

legend(2, 12, paste0(
  "Slope: ", signif(coef(res_lm)[2], 3),
  "\nBackground (mean): ", signif(res_bgpl[1], 3),
  "\nsd_bg: ", signif(res_sd_bg, 3),
  "\nPlateau: ", signif(res_bgpl[2], 3),
  "\ntop: ", signif(res_takeoff[[1]], 3),
  "\ntdp: ", signif(res_takeoff_reverse[[1]], 3)
), bty = "n")

mtext("A    Positive", cex = 1, side = 3, adj = 0, font = 2)


y_lim <- 2
plot(amp_data[, 1], amp_data[, 4], ylim = c(-0.1, y_lim), xlab = "Cycle", ylab = "RFU", type = "b", lwd = 2, pch = 19)
res_bgpl <- unlist(backgroundplateu(amp_data[, 4]))
abline(h = res_bgpl, col = c("green", "blue"))
abline(h = 0, col = "grey")

mtext("B    Negative", cex = 1, side = 3, adj = 0, font = 2)

curve_colors <- c(rainbow(ncol(boggy) - 1, alpha = .5))
matplot(boggy[, 1], boggy[, -1], type = "l", col = curve_colors, xlab = "Cycle", ylab = "RFU", lty = 1)
mtext("C    boggy dataset", cex = 1, side = 3, adj = 0, font = 2)

## ---- echo=FALSE, results='hide',message=FALSE, fig.cap="Commonly used methods for the analysis of quantification points. A) Linear plot of an amplification curve with a typical sigmoid shape. The Grey horizontal line is the threshold as determined by the \\textit{68-95-99.7 rule} from the fluoresce emission of cycle 1 to 10. The black horizontal line is the user defined threshold in the log-linear range of the amplification curve. The Ct is calculated from the intersection of the horizontal line and a quadratic polynomial fitted in to the amplification curve (see @roediger2015chippcr for details). B) The Amplification curve plot with a logarithmic ordinate visualizes the linear phase. C) Analysis of the amplification curve by fitting with a five parameter model (black line) (\\autoref{l5}). The red line is the first derivative of the amplification curve, with the maximum at 17.59 cycles. The maximum is used in selected system as Cq value and referred to as first derivative maximum (`cpD1`). The green line is the derivative of the amplification curve, with the maximum at 15.68 cycles a minimum approximately at 19.5 cycles. The maximum is used in selected system as Cq value and referred to as second derivative maximum (`cpD2`). The blue line is the amplification efficiency that is estimated from the trajactory of the exponential region. The `Eff` value 1.795 means that the amplification efficiency is approximately 89%. `cpDdiff` is the difference between the Cq values calculated from the first and the second derivative maximum ($cpDdiff = |cpD1 - cpD2|$) from the fitted model.\\label{figure_quntifcation_points}", fig.height=8, fig.width=8----
options(warn = -1)
library(qpcR)
library(chipPCR)
library(magrittr)

res_model <- pcrfit(testdat, cyc = 1, fluo = 2, model = l5)
res_takeoff <- takeoff(res_model, pval = 0.05, nsig = 3)

res_model_predict <- predict(res_model)

r_user <- 2.356

res_th.cyc <- th.cyc(testdat[, 1], testdat[, 2], r = r_user, linear = FALSE)


par(las = 0, oma = c(0, 0, 0, 0))

layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE), respect = TRUE)


plot(testdat[, 1], testdat[, 2], xlab = "Cycles", ylab = "Raw fluorescence")

abline(h = (mean(testdat[1:10, 2]) + 3 * sd(testdat[1:10, 2])), col = "grey")

abline(h = res_th.cyc[1, 2], col = "black")
text(28, r_user + 0.3, paste0("Threshold: ", r_user))
arrows(res_th.cyc[1, 1], res_th.cyc[1, 2], res_th.cyc[1, 1], 0, angle = 25, length = 0.1, lwd = 2)
mtext(paste0("A     ", "Ct = ", signif(res_th.cyc[1, 1], 4)), cex = 1.2, side = 3, adj = 0, font = 2)

plot(testdat[, 1], log(testdat[, 2]), xlab = "Cycles", ylab = "log(Raw fluorescence)")
abline(h = log(res_th.cyc[1, 2]), col = "black")
arrows(res_th.cyc[1, 1], log(res_th.cyc[1, 2]), res_th.cyc[1, 1], min(log(testdat[, 2]), na.rm = TRUE), angle = 25, length = 0.1, lwd = 2)
mtext(paste0("B     ", "Ct = ", signif(res_th.cyc[1, 1], 4)), cex = 1.2, side = 3, adj = 0, font = 2)

res_efficiency <- efficiency(res_model)
cpDdiff <- sqrt((res_efficiency$cpD1 - res_efficiency$cpD2)^2)

arrows(res_takeoff[[1]], res_takeoff[[2]], res_takeoff[[1]], -0.2, angle = 25, length = 0.1, lwd = 2)

abline(v = 19.5)

mtext(paste0("C   ",  "cpDdiff: ", cpDdiff), cex = 1.2, side = 3, adj = 0, font = 2)

## ---- echo=FALSE, fig.cap="Amplification curves of the `RAS002` dataset. All amplification curves of the `RAS002` dataset were manually classified (`negative`,`positive`). A) The negative amplification curves have no sigmoid curve progression. Two groups with different signal levels form the amplification curves. B) All positive amplification curves have a sigmoid curve shape and a similar ground signal. C) The density function of the RFU values from the first 15 PCR cycles shows a bimodal distribution. Based on these data, it is easy to divide them into two groups. Both groups' density functions appear to be symmetrical. D) The density function from the RFU values of the first 15 PCR cycles shows a monomodal distribution. It seems that the density function of the RFU values is left-skewed. \\label{amplification_curve_shapes}", fig.height=6----
options(warn = -1)
library(PCRedux)

index <- which(grepl("B.Globin", colnames(RAS002)))

data <- RAS002[, c(1, index)]

y_lim <- range(data[, -1])

curve_colors <- c(rainbow(ncol(data) - 1, alpha = .5))

upper_limit <- 15

cycles <- 2:upper_limit

par(mfrow = c(2, 2))

matplot(data[, 1], data[, which(RAS002_decisions[index - 1] == "n") + 1], type = "l", 
        col = curve_colors, xlab = "Cycle", ylab = "RFU", lty = 1, ylim = y_lim)
mtext("A    Negative", cex = 1.2, side = 3, adj = 0, font = 2)
abline(v = upper_limit)

matplot(data[, 1], data[, which(RAS002_decisions[index - 1] == "y") + 1], type = "l", 
        col = curve_colors, xlab = "Cycle", ylab = "RFU", lty = 1, ylim = y_lim)
mtext("B    Positive", cex = 1.2, side = 3, adj = 0, font = 2)
abline(v = upper_limit)


rfu_neg <- unlist(data[cycles, which(RAS002_decisions[index - 1] == "n") + 1])
plot(density(rfu_neg, width = 100), main = "", xlim = c(2300, 2950))
rug(rfu_neg)
mtext("C    Negative", cex = 1.2, side = 3, adj = 0, font = 2)


rfu_pos <- unlist(data[cycles, which(RAS002_decisions[index - 1] == "y") + 1])
plot(density(rfu_pos, width = 100), main = "", xlim = c(2300, 2950))
rug(rfu_pos)
mtext("D    Positive", cex = 1.2, side = 3, adj = 0, font = 2)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  # Install devtools, if not already installed.
#  install.packages("devtools")
#  
#  library(devtools)
#  install_github("devSJR/PCRedux")

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  # Select your local mirror
#  install.packages("PCRedux")

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  # The following command points to the help for download and install of packages
#  # from CRAN-like repositories or from local files.
#  ?install.packages()

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  library(PCRedux)
#  
#  context("qPCR2fdata")
#  
#  test_that("qPCR2fdata gives the correct dimensions and properties", {
#    library(qpcR)
#    res_fdata <- qPCR2fdata(testdat)
#    res_fdata_preprocess <- qPCR2fdata(testdat, preprocess = TRUE)
#  
#    expect_that(res_fdata, is_a("fdata"))
#    expect_that(length(res_fdata$rangeval) == 2 &&
#      res_fdata$rangeval[2] == 49, is_true())
#  
#    expect_that(res_fdata_preprocess, is_a("fdata"))
#    expect_that(length(res_fdata_preprocess$rangeval) == 2 &&
#      res_fdata_preprocess$rangeval[2] == 49, is_true())
#  })

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  library(RDML)
#  # Load the RDML package and use its functions to import the amplification curve
#  #  data
#  library(RDML)
#  filename <- system.file("RAS002.rdml", package = "PCRedux")
#  raw_data <- RDML$new(filename = filename)

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  # Export the RDML data from the PCRedux package as the objects RAS002 and RAS003.
#  library(RDML)
#  library(PCRedux)
#  library(magrittr)
#  suppressMessages(library(data.table))
#  
#  RAS002 <- data.frame(RDML$new(paste0(
#      path.package("PCRedux"),
#                                       "/", "RAS002.rdml"
#  ))$GetFData())
#  
#  # The obbject RAS002 can be stored in the working directory as CSV file with
#  # the name RAS002_amp.csv.
#  write.csv(RAS002, "RAS002_amp.csv", row.names = FALSE)

## ---- echo=TRUE, eval=TRUE, fig.cap="\\label{figure_curve_classification} Classification of amplification curves. The availability of classified amplification curves is an important prerequisite for the development of methods based on monitored learning. Amplification curves (n = 8858) from the `htPCR` dataset (\\texttt{qpcR} package, [@Ritz2008]) were classified in total eight time at different time points by a human eight times with the classes ambiguous (a), positive (y) or negative (n). The classification is subject to the subjectivity of the human expert, classified with the humanrater() function. Consequently, the amplification curves were selected randomly so that systematic errors in classification should be minimized. With this example, it becomes evident that even with the same dataset, different class assignments can occur. While in the first three rounds (A-C) only a few amplification curves were classified as negative. Their proportion is increased nearly tenfold (D-H) in subsequent classifications."----
# Suppress messages and load the packages for reading the data of the classified
# amplification curves.
options(warn = -1)
suppressMessages(library(data.table))
library(PCRedux)

# Load the decision_res_htPCR.csv dataset from a csv file.
filename <- system.file("decision_res_htPCR.csv", package = "PCRedux")
decision_htPCR <- fread(filename, data.table = FALSE)

#
par(mfrow = c(2, 4))
for (i in 2L:9) {
  data_tmp <- table(as.factor(decision_htPCR[, i]))

  barplot(data_tmp, col = adjustcolor("grey", alpha.f = 0.5),
      xlab = "Class", ylab = "Counts")
  text(c(0.7, 1.9, 3.1), rep(quantile(data_tmp, 0.25), 3), data_tmp, srt = 90)
  mtext(LETTERS[i - 1], cex = 1.2, side = 3, adj = 0, font = 2)
}

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  # Classify amplification curve data by correlation coefficients (r)
#  library(qpcR)
#  tReem(testdat[, 1:15], k = 3)

## ---- eval=TRUE, echo=FALSE, results='asis'------------------------------
options(warn = -1)
library(PCRedux)
library(xtable)
suppressMessages(library(data.table))
library(magrittr)

filename <- system.file("decision_res_htPCR.csv", package = "PCRedux")
decision_res_htPCR <- fread(filename, data.table = FALSE)

print(
  xtable(
    rbind(
      head(subset(decision_res_htPCR, conformity == FALSE), 25),
      head(subset(decision_res_htPCR, conformity == TRUE), 25)), 
         caption = "Results of the `htPCR` data set classification. 
         All amplification curves of the `htPCR` dataset were classified as 
         `negative`, `ambiguous` and `positive` by individuals in eight 
         analysis cycles (`test.result.1` $\\ldots$ `test.result.8`). If an 
         amplification curve has always been classified with the same class, 
         the last column (`conformity`) shows `TRUE`. As an example, the table 
         shows 25 amplification curves with consistent classes and 25 
         amplification curves with differing classes (`conformity = FALSE`).", 
         label = "tableheaddecision" ), include.rownames = FALSE, 
      comment = FALSE, caption.placement = "top", size = "\\tiny"
)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Use decision_modus() to go through each row of all classification done by
# a human.

dec <- lapply(1L:nrow(decision_res_htPCR), function(i) {
  decision_modus(decision_res_htPCR[i, 2:9])
}) %>% unlist()

names(dec) <- decision_res_htPCR[, 1]

# Show statistic of the decisions
summary(dec)

## ---- echo=FALSE, fig.cap="\\label{htPCR_nap}A) Comparison of amplification curves. Examples of a negative (black), ambiguous (red) and positive (green) amplification curve were selected from the `htPCR` dataset (\\texttt{qpcR} package, @Ritz2008). The negative amplification curve is not sigmoid and shows a strong positive trend. The ambiguous amplification curve approaches a sigmoid from, but shows a positive slope in the background (cycle 1 $\\rightarrow$ 5). The positive amplification curve is sigmoid. It begins in the background phase (cycle 1 $\\rightarrow$ 5) with a flat baseline, and shortly thereafter the exponential phase follows (cycle 5 $\\rightarrow$ 25) followed by a plateau phase (cycle 26 $\\rightarrow$ 35). B) Summary of frequencies of all classes of the `htPCR` record. negative, black; ambiguous, red; positive, green.", fig.height=5.5, fig.width=11----
options(warn = -1)
par(mfrow = c(1, 2))

library(qpcR)
matplot(
  htPCR[, 1], htPCR[, c(552, 512, 616)], xlab = "Cycle", ylab = "RFU",
  main = "htPCR dataset", type = "l", lty = 1, lwd = 2
)
legend("topleft", c(
  paste("negative ", colnames(htPCR)[552]),
  paste("ambiguos ", colnames(htPCR)[512]),
  paste("positive ", colnames(htPCR)[616])
), col = 1:3, pch = 19, bty = "n")
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)

# Plot the Frequencies of the decisions
barplot(
  table(dec), xlab = "Decision", ylab = "Frequency",
  main = "", col = c(2, 3, 1)
)
mtext("B    Classified by human", side = 3, adj = 0, font = 2, las = 0)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
library(PCRedux)
# Decisions for observation P01.W06
res_dec_P01.W06 <- decision_modus(decision_res_htPCR[
  which(decision_res_htPCR[["htPCR"]] == "P01.W06"),
  2L:9
], max_freq = FALSE)
print(res_dec_P01.W06)

## ---- echo=TRUE, fig.cap="\\label{visdat_pcrfit_plot}Application of ``visdat_pcrfit()`` for the visualization of the data structure after an analysis by ``pcrfit_single()``. The amplification curves (A01 = 1, A02 = 2, A04 = 3, B04 = 4) from the `C126EG685` dataset were analyzed with the ``pcrfit_single()`` function and then visualized with the ``visdat_pcrfit()`` function. For each observation, the classes (factor, integer, logical, numeric, NA) are presented. For the observations 2 and 4 the parameter `loglin_slope` could not be calculated (returned NA).", fig.height=8----
# Calculate curve features of an amplification curve dataset.
# Use the C126EG685 dataset from the chipPCR package and analyze the observations
# A01, A02, A04 and B05.

library(chipPCR)

res <- encu(C126EG685[, c(1,2,3,5,17)])

# Show all results in a plot. Note that the interactive parameter is set to
# FALSE.

visdat_pcrfit(res, type = "all", interactive = FALSE)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
# Calculate the Hausdorff distance of the amplification curves
# cluster the curves.
# Load additional packages for data and pipes.
options(warn = -1)
library(qpcR)
library(chipPCR)
suppressMessages(library(fda.usc))
library(magrittr)

# Convert the qPCR dataset to the fdata format
# Use unprocessed data from the testdat dataset
res_fdata <- qPCR2fdata(testdat)

# Extract column names and create rainbow color to label the data
columnames <- testdat[-1] %>% colnames()
data_colors <- rainbow(length(columnames), alpha = 0.5)

# Calculate the Hausdorff distance (fda.usc) package and plot the distances
# as clustered data.

res_fdata_hclust <- metric.hausdorff(res_fdata)
res_hclust <- hclust(as.dist(res_fdata_hclust))

# Cluster of the unprocessed amplification curves
res_cutree <- cutree(res_hclust, k = 2)
res_cutree <- factor(res_cutree)
levels(res_cutree) <- list(y = "1", n = "2")

## ---- echo=TRUE, fig.cap="\\label{qPCR2fdata}Shape based grouping of amplification curves. Grouping of amplification curves of the `testdat` dataset via Hausdorff distance. A) The amplification curves were converted with the qPCR2fdata() function. B) Subsequent they were processed by a cluster analysis using the Hausdorff distance. Faultless differentiation was achieved between negative amplification curves (n) and positive amplification curves (y).", fig.height=4, fig.width=8----
# Plot the converted qPCR data
par(mfrow = c(1, 2))
res_fdata %>% plot(
    ., xlab = "Cycle", ylab = "RFU", main = "", type = "l",
    lty = 1, lwd = 2, col = data_colors
)
legend(
    "topleft", paste0(as.character(columnames), ": ", res_cutree),
       pch = 19, col = data_colors, bty = "n", ncol = 2, cex = 0.7
)
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)

plot(res_hclust, main = "", xlab = "", sub="")
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)
rect(0.5, -3.5, 12.25, 0.5, border = "red")
text(7, 1, "negative", col = "red")
rect(12.5, -3.5, 24.5, 0.5, border = "green")
text(14, 1, "positive", col = "green", cex = 0.9)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
# Calculate slope and intercept on positive amplification curve data from the
# VideoScan 32 cavity real-time PCR device.
# Load additional packages for data and pipes.
options(warn = -1)
library(data.table)
library(fda.usc)
library(magrittr)

# Load the qPCR data from the HCU32_aggR.csv dataset
# Convert the qPCR dataset to the fdata format

filename <- system.file("HCU32_aggR.csv", package = "PCRedux")
data_32HCU <- fread(filename, data.table = FALSE)

res_fdata <- qPCR2fdata(data_32HCU)
# Extract column names and create rainbow color to label the data
columnames <- data_32HCU[-1] %>% colnames()
data_colors <- rainbow(length(columnames), alpha = 0.55)

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  # Load the qpcR package to calculate the Cq values by the second derivative
#  # maximum method.
#  options(warn = -1)
#  library(qpcR)
#  
#  res_Cq <- sapply(2L:ncol(data_32HCU), function(i) {
#      efficiency(pcrfit(data_32HCU, cyc = 1, fluo = i, model = l6))
#  })
#  
#  data.frame(
#      obs = colnames(data_32HCU)[-1],
#             Cq = unlist(res_Cq["cpD2", ]), eff = unlist(res_Cq["eff", ])
#  )
#  
#  #        Results
#  #
#  # obs    Cq      eff
#  # 1      A1 14.89 1.092963
#  # 2      B1 15.68 1.110480
#  # 3      C1 15.63 1.111474
#  # ...
#  # 30     F4 15.71 1.109634
#  # 31     G4 15.70 1.110373
#  # 32     H4 15.73 1.117827

## ---- echo=TRUE, fig.cap="\\label{HCU32}Clustering of amplification curves. The amplification curves from the 32HCU were processed with the ``qPCR2fdata()`` function and subsequent processed by a cluster analysis and Hausdorff distance analysis. A) Amplification curves were plotted from the raw data. B) Overall, the signal to noise ratios of the amplification curves were comparable between all cavities. C) The Cqs (Second Derivative Maximum) and the amplification efficiency (eff) were calculated with the ``efficiency(pcrfit())`` functions from the ``qpcR`` package. The median Cq is indicated as vertical line. Cqs larger or less than 0.1 of the Cq $\\tilde{x}$ are indicated with the labels of the corresponding observation. D) The clusters according to the Hausdorff distance show no specific pattern regarding the amplification curve signals. It appears that the observations D1, E1, F1, F3, G3 and H1 deviate most from the other amplification curves.", fig.height=8.5, fig.width=11----

library(fda.usc)
library(magrittr)

# To save computing time, the Cq values and amplification efficiencies were 
# calculated beforehand and transferred as a hard copy here.

calculated_Cqs <- c(
    14.89, 15.68, 15.63, 15.5, 15.54, 15.37, 15.78, 15.24, 15.94,
    15.88, 15.91, 15.77, 15.78, 15.74, 15.84, 15.78, 15.64, 15.61,
    15.66, 15.63, 15.77, 15.71, 15.7, 15.79, 15.8, 15.72, 15.7, 15.82,
    15.62, 15.71, 15.7, 15.73
)

calculated_effs <- c(
    1.09296326515231, 1.11047987547324, 1.11147389307153, 1.10308929700635,
    1.10012176315852, 1.09136717687619, 1.11871308210321, 1.08006168654712,
    1.09500422011318, 1.1078777171126, 1.11269436700649, 1.10628580163733,
    1.1082009954558, 1.11069683827291, 1.11074914659374, 1.10722949813473,
    1.10754282514113, 1.10098387264025, 1.1107026749644, 1.11599641663658,
    1.11388510347017, 1.11398547396991, 1.09410798249025, 1.12422338092929,
    1.11977386646464, 1.11212436173214, 1.12145338871426, 1.12180879952503,
    1.1080276005651, 1.10963449004393, 1.11037302758388, 1.11782689816295
)

# Plot the converted qPCR data
layout(matrix(c(1, 2, 3, 4, 4, 4), 2, 3, byrow = TRUE))
res_fdata %>% plot(
    ., xlab = "Cycle", ylab = "RFU", main = "HCU32_aggR", type = "l",
    lty = 1, lwd = 2, col = data_colors
)
legend(
    "topleft", as.character(columnames), pch = 19,
       col = data_colors, bty = "n", ncol = 4
)
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)

# Plot the background and plateau phase.

boxplot(
    data_32HCU[, -1] - apply(data_32HCU[, -1], 2, min),
        col = data_colors, las = 2, main = "Signal to noise ratio",
        xlab = "Sample", ylab = "RFU"
)
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)

# Plot the Cqs and the amplification efficiencies.
# Determine the median of the Cq values and label all Cqs, which a less 0.1 Cqs
# of the median or more then 0.1 Cqs of the median Cq.

plot(
    calculated_Cqs, calculated_effs, xlab = "Cq (SDM)",
     ylab = "eff", main = "Cq vs. Amplification Efficiency",
     type = "p", pch = 19, lty = 1, lwd = 2, col = data_colors
)

median_Cq <- median(calculated_Cqs)
abline(v = median_Cq)

text(median_Cq + 0.01, 1.085, expression(paste(tilde(x))))
labeled <- c(
    which(calculated_Cqs < median_Cq - 0.1),
             which(calculated_Cqs > median_Cq + 0.1)
)

text(
    calculated_Cqs[labeled], calculated_effs[labeled],
     as.character(columnames)[labeled]
)
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2)

# Calculate the Hausdorff distance using the fda.usc package and cluster the
# the distances.

res_fdata_hclust <- metric.hausdorff(res_fdata)
cluster <- hclust(as.dist(res_fdata_hclust))

# plot the distances as clustered data and label the leafs with the Cq values
# and colored dots.

plot(cluster, main = "Clusters of the amplification\n
curves as calculated by the Hausdorff distance", xlab = "", sub="")
mtext("D", cex = 1.2, side = 3, adj = 0, font = 2)

## ---- echo=TRUE, fig.cap="\\label{curve_fit_fail}Positive and negative Amplification curves from the `RAS002` dataset. A positive amplification curve (black) and a negative amplification curve (red) were selected from the `RAS002` dataset. The positive amplification curve has a baseline signal of approximately 2500 RFU and shows an unambiguous sigmoidal shape. The negative amplification curve has a baseline signal of approximately 4200 RFU, shows moderately positive slope and has no sigmoidal shape. B) A logistical function with seven parameters (`l7`) was adapted to the positive amplification curve. A Cq value of 25.95 was determined using the method of the maximum of the second derivative. The calculated Cq value appears to be correct. C) The negative amplification curve was also fitted with a seven parameter logistical function (`l7`). The second derivation method was used to determine the Cq value. Though a Cq value of 9.41 was calculated, it is clear that the model adaptation is not appropriate to calculate a trustworthy Cq value. If such a calculation would be done automatically, without human interaction, a false-positive result could be interpreted.", fig.height=7, fig.width=7----
# Load the qpcR package for the model fit.
suppressMessages(library(qpcR))
library(chipPCR)

# Select one positive and one negative amplification curve from the PCRedux 
# package.

amp_data <- RAS002[, c("cyc", "A01_gDNA.._unkn_B.Globin", "B07_gDNA.._unkn_HPRT1")]

# Arrange graphs in an matrix and set the plot parameters. An plot the positive
# and negative amplification curve.

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE), respect = TRUE)

matplot(amp_data[,1], amp_data[, -1], pch = 19, lty = 1, type = "l",
        xlab = "Cycle", ylab = "RFU", main = "")
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)

# Apply the the amptester function from the chipPCR package to the amplification 
# curve data and write the results to the main of the plots.

for (i in 2:3) {
    res.ampt <-  suppressMessages(amptester(amp_data[, i]))
    
    # Make a logical connection by two tests (shap.noisy, lrt.test and
    # tht.dec) of amptester to decide if an amplification reaction is
    # positive or negative.
    decision <- ifelse(!res.ampt@decisions[1] &&
    res.ampt@decisions[2] &&
    res.ampt@decisions[4],
    "positive", "negative"
    )
    # The amplification curves were fitted (l7 model) with pcrfit() function. The 
    # Cq was determined with the efficiency() function.
    
    fit <- pcrfit(data = amp_data, cyc = 1, fluo = i, model = l7)
    res <- efficiency(fit, plot = FALSE)
    plot(fit, pch = 19, lty = 1, type = "single", xlab = "Cycle", ylab = "RFU", 
         main = "", col = i - 1)
    abline(h = res[["fluo"]], v = res[["cpD2"]], col = c("grey", "red"))
    points(res[["cpD2"]], res[["fluo"]], pch = 19)
    
    mtext(paste0(LETTERS[i], "    Cq: ", res[["cpD2"]]), cex = 1.2, side = 3, 
          adj = 0, font = 2)
    
    legend(
        "topleft", paste0(colnames(amp_data)[i], "\nDecision: ", decision),
           bty = "n", cex = 1, col = "red"
    )
}

## ---- echo=TRUE----------------------------------------------------------
library(PCRedux)
str(pcrfit_single(RAS002[, 2]))

## ---- echo=FALSE, fig.cap="\\label{plot_models}Distribution of models of amplification curves. The `competimer`, `dil4reps94`, `guescini1`, `HCU32_aggR`, `karlen1`, `lc96_bACTXY`, `lievens1`, `RAS002`, `RAS003`, `reps384`, `rutledge`, `stepone_std`, `testdat`, `vermeulen1`, `VIMCFX96_60` datasets were analyzed with the ``encu()`` function. For each amplification curve, the optimal model was selected on the basis of the Akaike information criterion. A) Model functions of the raw amplification curve. B) Model functions of the rotated and flipped amplification curves. lNA, no model fitted. l4 \\ldots l7, model with four to seven parameters.", fig.height=3.5----
options(warn = -1)
library(PCRedux)

x <- data_sample$decision
y <- factor(data_sample[["qPCRmodel"]], levels = c("lNA", "l4", "l5", "l6", "l7"))

res_fw <- rbind(n = table(y[x == "n"]),
                y = table(y[x == "y"])
)

y <- factor(data_sample[["qPCRmodelRF"]], levels = c("lNA", "l4", "l5", "l6", "l7"))

res_rv <- rbind(n = table(y[x == "n"]),
                y = table(y[x == "y"])
)

colors <- c("grey", adjustcolor("red", alpha.f = 0.25))
par(mfrow = c(1,2))

barplot(res_fw / sum(res_fw) * 100, beside = TRUE, col = colors, xlab = "Model", ylab = "Percentage")
legend("top", legend = c("Negative", "Positive"), fill = colors, bty = "n")
mtext("A    Normal", cex = 1.2, side = 3, adj = 0, font = 2)

barplot(res_rv / sum(res_rv) * 100, beside = TRUE, col = colors, xlab = "Model", ylab = "Percentage")
legend("top", legend = c("Negative", "Positive"), fill = colors, bty = "n")
mtext("B    Rotated and Flipped", cex = 1.2, side = 3, adj = 0, font = 2)

## ---- echo=FALSE, fig.cap="\\label{plot_dat_Cq}Distribution of Cq values of positive and negative amplification curves. All Cq values were calculated from 3302 amplification curves after fitting the optimal multi-parametric models. The Cqs of positive amplification curves heaped up in the range between 10 and 35 PCR cycles. This differs from the distribution of negative amplification curves. The Cqs of negative amplification curves were calculate over the entire range. Note: The Cqs of the negative amplification curves are false positive. A) The maximum of the first derivative cpD1. B) The maximum of the second derivative (cpD2).", fig.height=5----
library(PCRedux)

data <- data_sample
feature <- c("cpD1", "cpD2")

par(mfrow = c(1,2))

x <- data$decision

for(i in 1L:length(feature)) {
    y <- data[, colnames(data) == feature[i]]
    y_density_neg <- density(y[x == "n"])
    y_density_pos <- density(y[x == "y"])
    res <- stats::wilcox.test(y ~ x)
    h <- max(na.omit(y))
    l <- min(na.omit(y))
    h_text <- rep(h * 0.976, 2)
    
    par(bg=NA)
    stripchart(y ~ x, vertical = TRUE, ylab = feature[i],
               method = "jitter", pch = 19, cex = 1, 
               col = adjustcolor("darkgrey", alpha.f = 0.65), 
               ylim = c(l * 0.95, h * 1.05))
    
    polygon(y_density_neg$y * 5 + 1.35, y_density_neg$x, col = adjustcolor("red", alpha.f = 0.25), border = NA)
    polygon(y_density_pos$y * 5 + 1.35, y_density_pos$x, col = adjustcolor("green", alpha.f = 0.25), border = NA)
    arrows(1,30,1.5,32, length = 0.1)
    arrows(2,5,1.5,7, length = 0.1)
    
    boxplot(y ~ x, outline = FALSE, add = TRUE, boxwex = 0.35)
    
    legend("topleft", paste0("P = ", signif(res[["p.value"]])), 
           cex = 1, bty = "n")
    
    mtext(paste0(LETTERS[i], "   ", feature[i]), cex = 1.2, side = 3, 
          adj = 0, font = 2, col = ifelse(signif(res[["p.value"]]) < 0.05, 
                                          "black", "red"))
}

## ---- echo=FALSE, results='hide',message=FALSE, fig.cap="Both the minimum (cpD2m) and the maximum (cpD2) of the second derivative were determined numerically using the ``diffQ2()`` function. In addition, the function returns the maximum of the first derivative (cpD1). The difference of cpD2 and cpD2m results in the `cpD2_range`. Large `cpD2_range` values indicate a low amplification efficiency or negative amplification reaction. `bg.start` provides an estimate for the end of the ground phase. The following formula is used for the calculation: $bg.start = cpD1 - f * (cpD2m - cpD2)$. The distance between cpD1 and cpD2 is multiplied by a factor. `bg.stop` provides an estimate for the end of the exponential phase. The following formula is used for the calculation: $bg.start = cpD1 + f * (cpD2m - cpD2)$. $f$ is a factor (default 0.6) (see manual of ``bg.max()`` for details).\\label{figure_cpD2_range}", fig.height=4----
library(chipPCR)
data <- cbind(RAS002[, 1], CPP(RAS002[, c(1,2)], method.norm = "minm")$y.norm)

# Invoke the inder function for the object data to interpolate 
# the derivatives of the simulated data as object res. The Nip 
# parameter was set to 5. This leads to smoother curves. res is
# an object of the class "der".
res <- inder(data[-c(1,2),], Nip = 5)

# Plot the object res and add descriptions to the elements.

par(las = 0, bty = "n", oma = c(.5,.5,.5,.5))

plot(data, xlab = "Cycle", ylab = "RFU", ylim = c(0,1),
     main = "", type = "b", pch = 20, lwd = 2)

# Add graphical elements for the derivatives and the calculated
# Cq values FDM, SDM and SDm.

lines(res[, "x"], res[, "d1y"], col = "blue", lwd = 2)
lines(res[, "x"], res[, "d2y"], col = "red", lwd = 2)

# Fetch the Cq values from res with the summary function
summ <- summary(res, print = FALSE)

abline(v = c(summ[[1]], summ[[2]], summ[[3]]), col = c("blue", "red", "red"), 
       lwd = 2)

text(summ[[2]] - 2, 0.5, paste("cpD2\n", round(summ["SDM"], 2)), 
     cex = 1.1, col = "red")
text(summ[[1]] + 2, 0.7, paste("cpD1\n", round(summ["FDM"], 2)), 
     cex = 1.1, col = "blue")
text(summ[[3]] + 2, 0.5, paste("cpD2m\n", round(summ["SDm"], 2)), 
     cex = 1.1, col = "red")

arrows(summ[[2]], 0.85, summ[[3]], 0.85, code = 3, length = 0.1)
text(summ[[1]] - 6, 0.85, paste("cpD2_range \n", 
                           round(abs(summ["SDM"] - summ["SDm"]), 2)), 
     cex = 1.1, col = "black")

legend("topleft", c("Amplification curve", "First derivative", "Second derivative"), 
       col = c(1,4,2), lty = c(2,1,1), bty = "n")

res <- data
background <- bg.max(res[-c(1,2), 1], res[-c(1,2), 2])


abline(v = background@bg.stop, col = "grey")
text(background@bg.stop, 0.3, "bg.stop", pos = 4, col = "grey", srt=90)
abline(v = background@amp.stop, col = "grey")
text(background@amp.stop, 0.3, "amp.stop", pos = 4, col = "grey", srt=90)

## ---- echo=FALSE, fig.cap="\\label{plot_dat_EffTop}Analysis of location features. Amplification curves from the datasets `stepone_std`, `RAS002`, `RAS003`, `lc96_bACTXY`, `C126EG595` and `dil4reps94` were analyzed with the ``encu()`` function. These datasets contain positive and negative amplification curves.  Furthermore, the meta dataset contains amplification curves that exhibit a hook effect or non-sigmoid shapes, for instance. All amplification curves are manually classified. Altogether 626 positive and 317 negative amplification curves were included in the analysis. A) `eff`, optimized PCR efficiency found within a sliding window. B) `sliwin`, PCR efficiency by the ‘window-of-linearity’ method. C) `cpDdiff`, difference between the Cq values calculated from the first and the second derivative maximum. D) `loglin_slope`, slope from the cycle at the second derivative maximum to the second derivative minimum. E) `cpD2_range`, absolute value of the difference between the minimum and the maximum of the second derivative maximum. F) `top`, takeoff point. G) `f.top`, fluorescence intensity at takeoff point. H) `tdp`,  takedown point. I) `f.tdp`, fluorescence intensity at takedown point. J) `bg.stop`, estimated end of the ground phase. K) `amp.stop`, estimated end of the exponential phase. L) `convInfo_iteratons`, number of iterations until convergence.", fig.height=11, fig.width=11----

library(PCRedux)

data <- data_sample[data_sample$dataset %in% c("stepone_std", "RAS002", "RAS003", "lc96_bACTXY", "C126EG595", "dil4reps94"), ]

feature <- c("eff", "sliwin", 
             "cpDdiff", "loglin_slope", "cpD2_range", "top", 
             "f.top", "tdp", "f.tdp", "bg.stop", 
             "amp.stop", "convInfo_iteratons")

par(mfrow = c(3,4))

x <- data$decision

for(i in 1L:length(feature)) {
    y <- data[, colnames(data) == feature[i]]
    res <- stats::wilcox.test(y ~ x)
    h <- max(na.omit(y))
    l <- min(na.omit(y))
    h_text <- rep(h * 0.976, 2)
    
    par(bg=NA)
    stripchart(y ~ x, vertical = TRUE, ylab = feature[i],
               method = "jitter", pch = 19, cex = 1, 
               col = adjustcolor("darkgrey", alpha.f = 0.65), 
               ylim = c(l * 0.95, h * 1.05))
    
    boxplot(y ~ x, outline = FALSE, add = TRUE, boxwex = 0.35)
    
    legend("topleft", paste0("P = ", signif(res[["p.value"]])), 
           cex = 1, bty = "n")
    
    mtext(paste0(LETTERS[i], "   ", feature[i]), cex = 1.2, side = 3, 
          adj = 0, font = 2, col = ifelse(signif(res[["p.value"]]) < 0.05, 
                                          "black", "red"))
}

## ---- eval=TRUE, echo=FALSE, fig.cap="\\label{loglin_slope}Concept of the `loglin_slope` feature. The algorithm determines the fluorescence values of the raw data at the approximate positions of the maximum of the first derivative, the minimum of the second derivative and the maximum of the second derivative, which are in the exponential phase of the amplification curve. A linear model is created from these parameter sets and the slope is determined. A) Positive amplification curves have a clearly positive slope. B) Negative amplification curves usually have a low, sometimes negative slope. The data were taken from the `RAS002` dataset.", fig.height=4----
# Load example data (observation F6.1) from the testdat dataset
options(warn = -1)
library(magrittr)
# Load MBmca package to calculate the minimum, and the maximum of the second
# derivative

library(MBmca)

par(mfrow = c(1, 2))

cyc <- 1

data <- RAS002[, c("cyc", "A01_gDNA.._unkn_B.Globin", "H10_ntc_ntc_B.Globin")]
col_names <- colnames(data)
data <- cbind(data[, 1], data[, 2]/quantile(data[, 2], .99), data[, 3]/quantile(data[, 3], .99))
colnames(data) <- col_names
# Calculate the minimum, and the maximum of the second
# derivative and assign it to the object res_diffQ2

rfu <- c(2,3)

for (i in rfu){
obs <- colnames(data)[i]
data_tmp <- data.frame(cyc = data[, cyc], rfu = data[, i])
res_diffQ2 <- suppressMessages(diffQ2(data_tmp, plot = FALSE, fct = min, inder = TRUE))
ROI <- round(c(res_diffQ2[[1]], res_diffQ2[[3]]))
# Build a linear model from der second derivative of res_diffQ2
res_loglin_lm <- lm(data_tmp[ROI, 2] ~ ROI)

data_tmp %>% plot(
    ., xlab = "Cycle", ylab = "normalized RFU", main = "", type = "l",
    lty = 1, lwd = 2, col = "black"
)
abline(res_loglin_lm, col = "red", lwd = 2)
abline(v = ROI, col = "grey", lwd = 2)
mtext(paste(LETTERS[i - 1], obs), side = 3, adj = 0, font = 1)
legend("topleft", paste("Slope: ", signif(
    coefficients(res_loglin_lm)["ROI"],
                                            4
)), bty = "n")
mtext("", side = 3, adj = 0, font = 2)
}

## ---- echo=FALSE, fig.cap="\\label{plot_sd_bg}Standard deviation in the ground phase of various qPCR devices. The `sd_bg` feature was used to determine if the standard deviation between the thermo-cyclers and between positive and negative amplification curves was different. The standard deviation was determined from the fluorescence values from the first cycle to the takeoff point. If the takeoff point could not be determined, the standard deviation from the first cycle to the eighth cycle was calculated. The Mann-Whitney test was used to compare the medians of the two populations (y, positive; n, negative). The differences were significant for A) LC_480 (Roche), B) ABI\\_Prism\\_7700 (ABI), C) LC1.0 (Roche), E) CFX96 (Bio-Rad) and F) LC96 (Roche). The difference was not significant for D) StepOne (Thermo Fisher).", fig.height=7----
library(PCRedux)
data <- data_sample[, c("device", "decision", "sd_bg")]

devices <- unique(data[data[["decision"]] == "n", "device"])

par(mfrow = c(3,3))
for(i in 1L:length(devices)) {
    data_tmp <- data[data[, 1] == devices[i], ]
    x <- data_tmp[["decision"]]
    y <- data_tmp[["sd_bg"]]
    pos_neg <- summary(x)
    c(paste0("y (", pos_neg[["y"]], ")"), paste0("n (", pos_neg[["n"]], ")"))
    res <- stats::wilcox.test(y ~ x)
    
    par(bg=NA)
    stripchart(y ~ x, vertical = TRUE, 
               xlab = "Decision",
               ylab = "Standard Deviation",
               method = "jitter", pch = 19, cex = 1, 
               col = adjustcolor("darkgrey", alpha.f = 0.65))
    
    boxplot(y ~ x, outline = FALSE, add = TRUE, boxwex = 0.35)
    legend("topleft", paste0("P = ", signif(res[["p.value"]])), 
           cex = 1, bty = "n")
    legend("bottomleft", c(paste0("y (", pos_neg[["y"]], ")"), 
                           paste0("n (", pos_neg[["n"]], ")")), 
           cex = 1, bty = "n")
    mtext(paste0(LETTERS[i], "   ", devices[i]), cex = 1.2, side = 3, 
          adj = 0, font = 2, col = ifelse(signif(res[["p.value"]]) < 0.05, 
                                          "black", "red"))
}

## ---- echo=FALSE, fig.cap="\\label{plot_bg_pt}Analysis of slope and ratio features. Amplification curves from the datasets `stepone_std`, `RAS002`, `RAS003`, `lc96_bACTXY`, `C126EG595` and `dil4reps94` were analyzed with the ``encu()`` function. These datasets contain positive and negative amplification curves.  Furthermore, the meta dataset contains amplification curves that exhibit a hook effect or non-sigmoid shapes, for instance. All amplification curves are manually classified. Altogether 626 positive and 317 negative amplification curves were included in the analysis. A) `b_slope`, B) `f_intercept`, C) `minRFU` is the minimum (1% qantile) of the amplification curve, D) `maxRFU` is the maximum (99% qantile) of the amplification curve, E) `init2` is the initial template fluorescence from an exponential model, F) `fluo` is the raw fluorescence value at the second derivative maximum, G) `slope_bg` is the slope calculated be the ``earlyreg()`` function, H) `intercept_bg` is the intercept calculated be the ``earlyreg()`` function, I) `sd_bg` is the standard deviation of the ground phase and J) `head2tail_ratio` is the between the RFU values of the head and the tail, normalized to the slope from the head to the tail.", fig.height=11, fig.width=11----
library(PCRedux)

data <- data_sample[data_sample$dataset %in% c("stepone_std", "RAS002", "RAS003", "lc96_bACTXY", "C126EG595", "dil4reps94"), ]

feature <- c("b_slope", "f_intercept", "minRFU", "maxRFU", 
             "init2", "fluo", "slope_bg", "intercept_bg", 
             "sd_bg", "head2tail_ratio")

par(mfrow = c(3,4))

x <- data$decision

for(i in 1L:length(feature)) {
    y <- data[, colnames(data) == feature[i]]
    res <- stats::wilcox.test(y ~ x)
    h <- max(na.omit(y))
    l <- min(na.omit(y)) 
    h_text <- rep(h * 0.976, 2)
    
    par(bg=NA)
    stripchart(y ~ x, vertical = TRUE, ylab = feature[i],
               method = "jitter", pch = 19, cex = 1, 
               col = adjustcolor("darkgrey", alpha.f = 0.65), 
               ylim = c(l * 0.95, h * 1.05))
    
    boxplot(y ~ x, outline = FALSE, add = TRUE, boxwex = 0.35)
    
    legend("topleft", paste0("P = ", signif(res[["p.value"]])), 
           cex = 1, bty = "n")
    
    mtext(paste0(LETTERS[i], "   ", feature[i]), cex = 1.2, side = 3, 
          adj = 0, font = 2, col = ifelse(signif(res[["p.value"]]) < 0.05, 
                                          "black", "red"))
}

## ---- echo=TRUE----------------------------------------------------------
options(warn = -1)
library(PCRedux)

data <- data_sample[data_sample$dataset %in% 
                c("batsch1",
                  "HCU32_aggR",
                  "lc96_bACTXY",
                  "RAS002", 
                  "RAS003", 
                  "stepone_std"), ]

n_positive <- sum(data[["decision"]] == "y")
n_negative <- sum(data[["decision"]] == "n")

dat <- data.frame(polyarea = data[, "polyarea"], 
                  decision = as.numeric(factor(data$decision, 
                                         levels = c("n", "y"), 
                                         label = c(0, 1))) - 1)

# Select randomly observations from 70% of the data for training.
# n_train is the number of observations used for training.

n_train <- round(nrow(data) * 0.7)

# index_test is the index of observations to be selected for the training
index_test <- sample(1L:nrow(dat), size = n_train)

# index_test is the index of observations to be selected for the testing
index_training <- which(!(1L:nrow(dat) %in% index_test))

# train_data contains the data used for training

train_data <- dat[index_test, ]

# test_data contains the data used for training

test_data <- dat[index_training, ]

# Fit the binomial logistic regression model

model_glm <- glm(decision ~ polyarea, family=binomial(link='logit'), 
                 data = train_data)

predictions <- ifelse(predict(model_glm, 
                                 newdata = test_data, type = 'response') > 0.5,
                         1, 0)

res_performeR <- performeR(predictions, test_data[["decision"]])[, c(1:10, 12)]


## ---- echo=TRUE----------------------------------------------------------
summary(model_glm)

## ---- echo=TRUE, fig.cap="\\label{plot_Logistic_Regression}Binomial logistic regression for the `polyarea` feature. A) binomial logistic regression model for the response variable $Y$ (decision) is categorical and must be converted into a numerical value. This regression calculation makes it possible to estimate the probability of a categorical response using predictor variables $X$. In this case, the predictor variable is `polyarea`. Gray dots are the value values used for training. Red dots are the values used for testing. The regression curve of the binomial logistic regression is shown in blue. At 0.5, the gray horizontal line marks the threshold value of probability used to determine whether an amplification curve is negative or positive. B) The measure were determined with the \\textit{performance()} function from the \\texttt{PCRedux} package. Sensitivity, TPR; Specificity, SPC; Precision, PPV; Negative predictive value, NPV; Fall-out, FPR; False negative rate, FNR; False discovery rate, FDR; Accuracy, ACC; F1 score, F1; Matthews correlation coefficient, MCC, Cohen's kappa (binary classification), kappa ($\\kappa$)."----
options(warn = -1)
library(PCRedux)

par(mfrow = c(1,2))

# Plot train_data (grey points) and the predicted model (blue)

plot(train_data$polyarea, train_data$decision, pch = 19, 
     xlab = "polyarea", ylab = "Probability", 
     col = adjustcolor("grey", alpha.f = 0.9), cex = 1.5)
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)
abline(h = 0.5, col = "grey")

curve(predict(model_glm, data.frame(polyarea = x), type = "resp"), 
      add = TRUE, col = "blue")

# Plot test_data (red)

points(test_data$polyarea, test_data$decision, pch = 19,
       col = adjustcolor("red", alpha.f = 0.3))
legend("right", paste("Positive: ", n_positive, 
                      "\nNegative: ", n_negative), bty = "n")


# Plot the sensitivity, specificity and other measures to describe the prediction.

position_bp <- barplot(as.matrix(res_performeR), yaxt = "n", 
                       ylab = "Probability", main = "", las = 2, 
                       col = adjustcolor("grey", alpha.f = 0.5))

par(srt = 90)
text(position_bp, rep(0.8, length(res_performeR)), 
     paste(signif(res_performeR, 2)*100, "%"), cex = 0.6)
axis(2, at = c(0, 1), labels = c("0", "1"), las = 2)
abline(h = 0.85, col = "grey")

mtext("B", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)

## ---- eval=TRUE, echo=FALSE----------------------------------------------
data(testdat, package = "qpcR")
x <- testdat[, 1]
y.pos <- testdat[, 2]
y.neg <- testdat[, 4]
nh <- trunc(length(x) * 0.20)
nt <- trunc(length(x) * 0.15)

y.pos.head <- head(y.pos, n = nh)
y.neg.head <- head(y.neg, n = nh)
y.pos.tail <- tail(y.pos, n = nt)
y.neg.tail <- tail(y.neg, n = nt)

lb.pos <- median(y.pos.head) + 2 * mad(y.pos.head)
ub.pos <- median(y.pos.tail) - 2 * mad(y.pos.tail)

lb.neg <- median(y.neg.head) + 2 * mad(y.neg.head)
ub.neg <- median(y.neg.tail) - 2 * mad(y.neg.tail)

res.shapiro.pos <- shapiro.test(y.pos)
res.shapiro.neg <- shapiro.test(y.neg)

res.wt.pos <- wilcox.test(head(y.pos, n = nh), tail(y.pos, n = nt), alternative = "less")
res.wt.neg <- wilcox.test(head(y.neg, n = nh), tail(y.neg, n = nt), alternative = "less")

###
RGt <- function(y) {
  ws <- ceiling((15 * length(y)) / 100)
  if (ws < 5) {
    ws <- 5
  }
  if (ws > 15) {
    ws <- 15
  }
  y.tmp <- na.omit(y[-c(1:5)])
  x <- 1:length(y.tmp)
  suppressWarnings(
    res.reg <- sapply(1L:(length(y.tmp)), function(i) {
      round(summary(lm(y.tmp[i:c(i + ws)] ~ x[i:c(i + ws)]))[["r.squared"]], 4)
    })
  )

  # Binarize R^2 values. Everything larger than 0.8 is positive
  res.LRt <- res.reg
  # Define the limits for the R^2 test
  res.LRt[res.LRt < 0.8] <- 0
  res.LRt[res.LRt >= 0.8] <- 1
  # Seek for a sequence of at least six positive values (R^2 >= 0.8)
  # The first five measuring point of the amplification curve are skipped
  # because most technologies and probe technologies tend to overshot
  # in the start (background) region.
  res.out <- sapply(5L:(length(res.LRt) - 6), function(i) {
    ifelse(sum(res.LRt[i:(i + 4)]) == 5, TRUE, FALSE)
  })
  out <- cbind(1L:(length(y.tmp)), res.reg)
  # res.out
}

## ---- echo = FALSE, fig.cap="The positive amplification curve F1.1 (`testdat` dataset) was analyzed with algorithms of the ``amptester()`` function. A) The Threshold test (THt) is based on the Wilcoxon rank sum test. The test compares 20% of the head to 15% of the tail region. A significant difference ($p-value = 0.000512$) between the two regions was found for the amplfication curve F1.1. This is indicative of a positive amplification reaction. B) Quantile-Quantile plot (Q-Q plot) of the the amplification curve. A Q-Q plot is a probability plot for a graphical comparison of two probability distributions by plotting their quantiles against each other. In this study the probability distribution of the amplification curve is compared to a theoretical normal distribution. The orange line is the theoretical normal quantile-quantile plot which passes through the probabilities of the first and third quartiles. The Shapiro-Wilk test (SHt) of normality checks whether the underlying population of a sample (amplification curve) is significantly ($\\alpha \\leq 5e^{-4}$) normal distributed. Since the p-value is $7.09 e^{-9}$ the null hypothesis can be rejected. C) The Linear Regression test (LRt). This test determines the coefficient of determination ($R^{2}$) by an ordinary least squares linear (OLS) regression. Usually the non-linear part of an amplification curve has an $R^{2}$ smaller than 0.8.\\label{statistical_methods_positive}", fig.height=7----

y_lim <- range(testdat[, 2]) * c(1, 1.3)

layout(matrix(c(1, 1, 2, 2, 1, 1, 3, 3), 2, 4, byrow = TRUE), respect = TRUE)

plot(
  x, y.pos, ylim = y_lim, xlab = "Cycle",
  ylab = "RFU", main = "", type = "b", pch = 19
)
mtext("A    Positive amplification", side = 3, adj = 0, font = 2)
abline(v = c(nh, length(x) - nt), lty = 3)

abline(h = lb.pos, lty = 2, col = "red")
text(35, 1.5, "Noise\nmedian + 2 * MAD", col = "red", cex = 1)

abline(h = ub.pos, lty = 2, col = "green")
text(35, 10, "Signal\nmedian - 2 * MAD", col = "green", cex = 1)

arrows(4.5, 12.5, 42.5, 12.5, length = 0.1, angle = 90, code = 3)
text(25, 14.5, paste(
  "W=", signif(res.wt.pos$statistic, 3), "\np-value = ",
  signif(res.wt.pos$p.value, 3)
))
# text(25, 5, paste("Fold change: \n", round(ub.pos / lb.pos, 2)))

qqnorm(y.pos, pch = 19, main = "")
legend("bottomright", paste(
  "SHt, W = ", signif(res.shapiro.pos$statistic, 3),
  "\np-value = ", signif(res.shapiro.pos$p.value, 3)
), bty = "n")
mtext("B", side = 3, adj = 0, font = 2)
qqline(y.pos, col = "orange", lwd = 2)

plot(
  RGt(y.pos), xlab = "Cycle", ylab = expression(R ^ 2), main = "", pch = 19,
  type = "b"
)
mtext("C    LRt", side = 3, adj = 0, font = 2)
abline(h = 0.8, col = "black", lty = 2)

## ---- echo = FALSE, fig.cap="The negative amplification curve F1.3 (`testdat` dataset) was analyzed with algorithms of the ``amptester()`` function. A) The Threshold test (THt) is based on the Wilcoxon rank sum test. The test compares 20% of the head to 15% of the tail region. No significant difference between the two regions was found for the amplfication curve F1.3. Since the p-value is $0.621$ the null hypothesis cannot be rejected. This is indicative of a negative amplification reaction. B) Quantile-Quantile plot (Q-Q plot) of the the amplification curve. A Q-Q plot is a probability plot for a graphical comparison of two probability distributions by plotting their quantiles against each other. In this study the probability distribution of the amplification curve is compared to a theoretical normal distribution. The orange line is the theoretical normal quantile-quantile plot which passes through the probabilities of the first and third quartiles. The Shapiro-Wilk test (SHt) of normality checks whether the underlying population of a sample (amplification curve) is significantly ($\\alpha \\leq 5e^{-4}$) normal distributed. Since the p-value is $0.895$ the null hypothesis cannot be rejected. C) The Linear Regression test (LRt). This test determines the coefficient of determination ($R^{2}$) by an ordinary least squares linear (OLS) regression. Usually the non-linear part of an amplification curve has an $R^{2}$ smaller than 0.8.\\label{statistical_methods_negative}", fig.height=7----
layout(matrix(c(1, 1, 2, 2, 1, 1, 3, 3), 2, 4, byrow = TRUE), respect = TRUE)
y_lim <- range(testdat[, 4]) * c(1.5, 2.5)

plot(
  x, y.neg, ylim = y_lim, xlab = "Cycle",
  ylab = "RFU", main = "", type = "b", pch = 19
)
mtext("B    Negative amplification", side = 3, adj = 0, font = 2)
abline(v = c(nh, length(x) - nt), lty = 3)

abline(h = lb.neg, lty = 2, col = "red")
text(10, 0.04, "Noise\nmedian + 2 * MAD", col = "red", cex = 1)

abline(h = ub.neg, lty = 2, col = "green")
text(40, -0.02, "Signal\nmedian - 2 * MAD", col = "green", cex = 1)

arrows(4.5, 0.04, 42.5, 0.04, length = 0.1, angle = 90, code = 3)
text(10, 0.06, paste(
  "W = ", res.wt.neg$statistic, "\np-value = ",
  signif(res.wt.neg$p.value, 3)
))
# legend(
#   "topright", paste("Fold change: \n", round(ub.neg / lb.neg, 2)),
#   bty = "n"
# )


qqnorm(y.neg, pch = 19, main = "")
legend(
  "bottomright", paste(
    "SHt, W = ", signif(res.shapiro.neg$statistic, 3),
    "\np-value = ", signif(res.shapiro.neg$p.value, 3)
  ),
  bty = "n"
)
mtext("B", side = 3, adj = 0, font = 2)
qqline(y.neg, col = "orange", lwd = 2)

plot(
  RGt(y.neg), xlab = "Cycle", ylab = expression(R ^ 2), main = "",
  pch = 19, type = "b"
)
mtext("C    LRt", side = 3, adj = 0, font = 2)
abline(h = 0.8, col = "black", lty = 2)

## ---- echo=FALSE, eval=FALSE, fig.cap="\\label{figure_autocorrelation_tau} Effect of tau.", fig.height=5----
#  # library(PCRedux)
#  #
#  # amp_data <- RAS002
#  #
#  # tau <- seq(1, 25, 1)
#  #
#  # par(mfrow = c(2, 2))
#  #
#  # plot(amp_data[, 1], amp_data[, 10], xlab = "Cycle", ylab = "RFU", pch = 19)
#  # mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)
#  #
#  # ac_pos <- sapply(1L:length(tau), function(i) {
#  #   autocorrelation_test(amp_data[, 2], n = tau[i], ns_2_numeric = TRUE)
#  # })
#  #
#  # plot(tau, ac_pos, ylim = c(0, 1), pch = 19)
#  # points(
#  #   tau[ac_pos == min(ac_pos)], ac_pos[ac_pos == min(ac_pos)], pch = 1, cex = 2,
#  #   col = "red"
#  # )
#  # legend(
#  #   "bottomleft", paste(
#  #     "tau: ", tau[ac_pos == min(ac_pos)],
#  #     "\nAC: ", signif(ac_pos[ac_pos == min(ac_pos)])
#  #   ),
#  #   bty = "n"
#  # )
#  # mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)
#  #
#  # plot(amp_data[, 1], amp_data[, 15], xlab = "Cycle", ylab = "RFU", pch = 19)
#  # mtext("C", cex = 1.2, side = 3, adj = 0, font = 2)
#  #
#  # ac_neg <- sapply(1L:length(tau), function(i) {
#  #   autocorrelation_test(amp_data[, 15], n = tau[i])
#  # })
#  #
#  # plot(tau, ac_neg, ylim = c(0, 1), pch = 19)
#  #
#  # mtext("D", cex = 1.2, side = 3, adj = 0, font = 2)

## ---- echo=TRUE, fig.cap="\\label{autocorrelation}Autocorrelation analysis for amplification curves of the `RAS002` dataset (\\texttt{PCRedux} package). A) Plot of all amplification curves of the `RAS002` dataset. B) Density plot of B) Positive curves and negative curves as determined by the `autocorrelation_test()` and a human expert. C) Performance analysis by the ``performeR()`` function (see \\autoref{section_performeR} for details).", fig.height=5.5, fig.width=11----
# Test for autocorrelation in amplification curve data
# Load the libraries magrittr for pipes and the amplification curve the data
# The amplification curve data from the `RAS002` dataset was used.
# The data.table package was used for fast import of the csv data
options(warn = -1)
library(magrittr)
library(PCRedux)
suppressMessages(library(data.table))

data <- RAS002

# Test for autocorrelation in the RAS002 dataset

res_ac <- sapply(2:ncol(data), function(i) {
  autocorrelation_test(data[, i], ns_2_numeric = TRUE)
})

# Curves classified by a human after analysis of the overview. 1 = positive,
# 0 = negative

human_classification <- fread(system.file(
  "decision_res_RAS002.csv",
  package = "PCRedux"
))

head(human_classification)


decs <- sapply(1L:nrow(human_classification), function(i) {
  res <- decision_modus(human_classification[i, 2L:(ncol(human_classification) - 1)])
  if (length(res) > 1) res[[1]] <- "n"
  res[[1]]
}) %>% unlist()


# Plot curve data as overview
# Names of the observations

layout(matrix(c(1, 2, 3, 1, 4, 4), 2, 3, byrow = TRUE))
matplot(
  data[, 1], data[, -1], xlab = "Cycle", ylab = "RFU",
  main = "", type = "l", lty = 1,
  col = decs, lwd = 2
)
legend("topleft", c("positive", "negative"), pch = 19, col = c(1, 2), bty = "n")
mtext("A    RAS002 dataset", cex = 1.2, side = 3, adj = 0, font = 2)


# Convert the n.s. (not significant) in 0 and others to 1.
# Combine the results of the aromatic autocorrelation_test as variable "ac",
# the human classified values as variable "hc" in a new data frame (res_ac_hc).

cutoff <- 0.8

res_ac_hc <- data.frame(
  ac = ifelse(res_ac > cutoff, 1, 0),
  hc = as.numeric(as.factor(decs)) - 1
) %>% as.matrix()
res_performeR <- performeR(res_ac_hc[, "ac"], res_ac_hc[, "hc"])


plot(density(res_ac), ylab = "Autocorrelation", main = "")
rug(res_ac)

abline(v = cutoff)
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)


cdplot(
  as.factor(decs) ~ res_ac, xlab = "Autocorrelation",
  ylab = "Decision"
)
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)

barplot(
  as.matrix(res_performeR[, c(1:10, 12)]), yaxt = "n", ylab = "",
  main = "Performance of autocorrelation_test",
  col = adjustcolor("grey", alpha.f = 0.5)
)

axis(2, at = c(0, 1), labels = c("0", "1"), las = 2)
mtext("D", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)

## ---- echo = TRUE, fig.cap = "\\label{earlyreg_slopes}Analysis of the ground phase with the ``earlyreg()`` function. A) The amplification curves show different slopes and intercepts in the early ground phase (ROI: cycle 1 to 5) of the qPCR. Amplification curves (n = 192) from the `RAS002` dataset were used. B) Both the slope and the intercept were used for a cluster analysis (k-means, Hartigan-Wong algorithm, number of centers \\textit{k = 5}). The amplification curves were separated into three clusters dependent on their slope and intercept (colored in red, green, cyan, balck). ", fig.height=7----
options(warn = -1)
library(PCRedux)

data <- RAS002

well <- substr(colnames(data)[-1], 1, 10)


# Normalize each amplification curve to their 0.99 percentile and use the
# earlyreg function to determine the slope and intercept of the first
# 5 cycles

res_earlyreg <- do.call(rbind, lapply(2L:ncol(data), function(i) {
  earlyreg(x = data[, 1], y = data[, i], range = 5, normalize = FALSE)
}))

# Label the observation with their original names
rownames(res_earlyreg) <- colnames(data)[2:ncol(data)]


cl <- kmeans(res_earlyreg, 5)

rownames(res_earlyreg) <- well

par(fig = c(0,1,0,1), las = 0, bty = "o", oma = c(0, 0, 0, 0))
matplot(
  data[, 1], data[, -1], pch = 19, lty = 1, type = "l",
  xlab = "Cycle", ylab = "RFU", main = "", col = cl[["cluster"]]
)
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)
abline(v = c(1,5))
rect(20.5,3500,45,4700, col = "white", border = NA)
text(3, 3250, "ROI")

par(fig = c(0.525, 0.99, 0.5, 0.95), new = TRUE) 
plot(res_earlyreg, col = cl[["cluster"]], pch = 19)
mtext("B    k-means, k = 5", cex = 1.2, side = 3, adj = 0, font = 2)

## ---- echo=TRUE, fig.cap="\\label{figure_head2tailratio}Calculation of the ratio between the head and the tail of a quantitative PCR amplification curve. A) Plot of quantile normalized amplification curves from the `RAS002` dataset. ROIs of the head and and tail are highlighted by circles. The ranges for performing Robust Linear Regression are automatically selected using the 25% and 75% quantiles. Therefore not all data points are used in the regression model. The straight line is the regression line from the robust linear model. The slopes of the positive and negative amplification curves differ. B) Boxplot for the comparison of the head to tailratio. Positive amplification curves have a lower ratio than negative curves. The difference between the classes is significant. ", fig.height=5.5, fig.width=11----
options(warn = -1)
library(PCRedux)

# Load the RAS002 dataset and assign it to the object data

data <- RAS002
data_decisions <- RAS002_decisions

# Calculate the head2tailratio of all amplification curves

res_head2tailratio <- lapply(2L:ncol(data), function(i) {
  head2tailratio(
    y = data[, i], normalize = TRUE, slope_normalizer = TRUE,
    verbose = TRUE
  )
})

# Fetch all values of the head2tailratio analysis for a later comparison
# by a boxplot.

res <- sapply(1L:length(res_head2tailratio), function(i)
  res_head2tailratio[[i]]$head_tail_ratio)

data_normalized <- cbind(
  data[, 1],
  sapply(2L:ncol(data), function(i) {
    data[, i] / quantile(data[, i], 0.99)
  })
)

# Assign color to the positive and negative decisions
colors <- as.character(factor(
  data_decisions, levels = c("y", "n"),
  labels = c(
      adjustcolor("black", alpha.f = 0.25), adjustcolor("red", alpha.f = 0.25))
))

res_wilcox.test <- stats::wilcox.test(res ~ data_decisions)

h <- max(na.omit(res))
h_text <- rep(h * 0.976, 2)

# Plot the results of the analysis

par(mfrow = c(1, 2), las = 0, bty = "o", oma = c(0, 0, 0, 0))

matplot(
  data_normalized[, 1], data_normalized[, -1],
  xlab = "Cycle", ylab = "normalized RFU", main = "RAS002 dataset",
  type = "l", lty = 1, lwd = 2, col = colors
)
for (i in 1L:(ncol(data_normalized) - 1)) {
  points(
    res_head2tailratio[[i]]$x_roi, res_head2tailratio[[i]]$y_roi,
    col = colors[i], pch = 19, cex = 1.5
  )
  abline(res_head2tailratio[[i]]$fit, col = colors[i], lwd = 2)
}
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)

# Boxplot of the head2tail ratios of the positive and negative
# amplification curves.

boxplot(res ~ data_decisions, col = unique(colors), ylab = "Head to Tail Ratio")

lines(c(1, 2), rep(h * 0.945, 2))
text(1.5, h_text, paste0("P = ", signif(res_wilcox.test[["p.value"]])), 
     cex = 1)

mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Calculate slope and intercept on noise (negative) amplification curve data
# for the last eight cycles.
options(warn = -1)
library(qpcR)
library(magrittr)

res_hook <- sapply(2:ncol(boggy), function(i) {
  hookreg(x = boggy[, 1], y = boggy[, i])
}) %>%
  t() %>%
  data.frame(obs = colnames(boggy)[-1], .)

## ---- eval=TRUE, echo=FALSE, results='asis'------------------------------
library(xtable)
print(
  xtable(
    res_hook, caption = "Screening results for the analysis with the 
             hookreg algorithm. Samples withe a value of 1 in the hook column 
             had all a hook effect like curvature. The observations F4.1, F4.2, F5.1, 
             F5.2, F6.1 and F6.2 miss entries because the hoogreg algorithm 
             could not fit a linear model. This is an expected behavior, since 
             these amplification curves did not have a hook effect like 
             curvature.",
    label = "tablehookreg"
  ), include.rownames = FALSE, comment = FALSE, caption.placement = "top",
  size = "\\tiny"
)

## ---- echo=TRUE, fig.cap="\\label{plot_hookreg}Detection of the hook effect in amplification curves. Amplification curves of the `boggy` dataset (\\texttt{qpcR}) were analyzed using the ``hookreg()`` function. The hook effect is characterized by a negative slope in the supposed plateau phase. Samples F1.1, F1.2, F2.1 and F2.2 (red) show a statistically significant negative slope. In the red rectangle the area of the hook effect is highlighted exemplarily for the sample F1.1. No statistically significant negative slope (no hook effect) could be observed for the remaining amplification curves (black).", fig.height=5, fig.width=7----
matplot(
  x = boggy[, 1], y = boggy[, -1], xlab = "Cycle", ylab = "RFU",
  main = "", type = "l", lty = 1, lwd = 2, col = res_hook$hook + 1
)

res_hook$hook.start[res_hook$hook.start == 0] <- "no hook"

legend(
  "topleft", paste(as.character(res_hook$obs), "Hook start cycle:", res_hook$hook.start), pch = 19,
  col = res_hook$hook + 1, bty = "n"
)
rect(26,2,40,2.3, border ="red")

## ---- echo=TRUE, eval = TRUE, fig.cap="\\label{plot_mblrr}Robust local regression to analyze amplification curves. The amplification curves were arbitrarily selected from the `RAS002` dataset. Not the differences in slopes and intercepts (red and green lines). The ``mblrr()`` function is presumably useful for datasets which are accompanied by noise and artifacts.m, slop; n, intercept.", fig.height=6----
options(warn = -1)
library(PCRedux)

# Select four amplification curves from the RAS002 dataset

data <- RAS002[, c(1, 2, 3, 4, 5)]


par(mfrow = c(2, 2))

for (i in 2L:ncol(data)) {
  x <- data[, 1]
  y_tmp <- data[, i] / quantile(data[, i], 0.99)
  res_q25 <- y_tmp < quantile(y_tmp, 0.25)
  res_q75 <- y_tmp > quantile(y_tmp, 0.75)
  res_q25_lm <- try(
    suppressWarnings(lmrob(y_tmp[res_q25] ~ x[res_q25])),
    silent = TRUE
  )
  res_q75_lm <- try(
    suppressWarnings(lmrob(y_tmp[res_q75] ~ x[res_q75])),
    silent = TRUE
  )

  plot(x, y_tmp, xlab = "Cycle", ylab = "RFU (normalized)",
    main = "", type = "b", pch = 19)
  
  mtext(paste0(LETTERS[i], "   ", colnames(data)[i]), cex = 1, side = 3, 
        adj = 0, font = 2)
  abline(res_q25_lm, col = "red")
  points(x[res_q25], y_tmp[res_q25], cex = 2.5, col = "red")
  abline(res_q75_lm, col = "green")
  points(x[res_q75], y_tmp[res_q75], cex = 2.5, col = "green")
}

## ---- eval=TRUE, echo=TRUE, results='asis'-------------------------------
# Load the xtable library for an appealing table output
library(xtable)

# Analyze the data via the mblrr() function

res_mblrr <- do.call(cbind, lapply(2L:ncol(data), function(i) {
  suppressMessages(mblrr(
    x = data[, 1], y = data[, i],
    normalize = TRUE
  )) %>% data.frame()
}))
colnames(res_mblrr) <- colnames(data)[-1]

# Transform the data for a tabular output and assign the results to the object
# output_res_mblrr.

output_res_mblrr <- res_mblrr %>% t()

# The output variable names of the mblrr() function are rather long. For better
# readability the variable names were changed to "nBG" (intercept of head region),
# "mBG" (slope of head region), "rBG" (Pearson correlation of head region),
# "nTP" (intercept of tail region), "mTP" (slope of tail region), "rBG" (Pearson
# correlation of tail region)

colnames(output_res_mblrr) <- c(
  "nBG", "mBG", "rBG",
  "nTP", "mTP", "rTP"
)

print(xtable(
  output_res_mblrr, caption = "mblrr() text intro. nBG, intercept of 
             head region; mBG, slope of head region; rBG, Pearson 
             correlation of head region; nTP, intercept of tail region; mTP, 
             slope of tail region; rBG, Pearson correlation of tail region",
  label = "tablemblrrintroduction"
), comment = FALSE, caption.placement = "top")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
# Load the xtable library for an appealing table output
options(warn = -1)
suppressMessages(library(FFTrees))
library(PCRedux)

# The RAS002 amplification curves were analyzed with the mblrr() function 
# to save computing time and the.results of this analysis are stored in the 
# `data_sample` dataset.

data <- data_sample[data_sample$dataset == "RAS002", c("mblrr_intercept_bg", 
                                                       "mblrr_slope_bg", 
                                                       "mblrr_cor_bg", 
                                                       "mblrr_intercept_pt", 
                                                       "mblrr_slope_pt", 
                                                       "mblrr_cor_pt")]

# The output variable names of the mblrr() function are rather long. For better
# readability the variable names were changed to "nBG" (intercept of head
# region), "mBG" (slope of head region), "rBG" (Pearson correlation of head
# region), "nTP" (intercept of tail region), "mTP" (slope of tail region),
# "rBG" (Pearson correlation of tail region).


res_mblrr <- data.frame(
    class = as.numeric(as.character(factor(RAS002_decisions, 
                                           levels = c("y", "n"), 
                                           label = c(1, 0)))),
  data
)

colnames(res_mblrr) <- c("class", "nBG", "mBG", "rBG", "nTP", "mTP", "rTP")

res_mblrr.fft <- suppressMessages(
            FFTrees(formula = class ~., data = res_mblrr)
            )

## ---- echo=FALSE, fig.cap="\\label{plot_FFTrees}Visualization of FFTrees of a ``mblrr()`` function analysis. **Top row** `Data`) Overview of the dataset, with displaying the total number of observations (N = 192) and percentage of positive (22%) and negative (78%) amplification curves. **Middle row** `FFT #1 (of 6)`) Decision Tree with the number of observations classified at each level of the tree. For the analysis, six features (nBG, intercept of head region; mBG, slope of head region; rBG, Pearson correlation of head region; nTP, intercept of tail region; mTP, slope of tail region; rBG, Pearson correlation of tail region) have been used for the analysis. After two tree levels (nBG, nTP) already the decision tree is created. All positive amplification curves (N = 40) are correctly classified. Two observations are classified as false-negative in the negative amplification curves. **Lower row** `Performance`)  The ``FFTrees()`` function determines several performance statistics. For the training data, there is a classification table on the left side showing the relationship between tree `decision` and the `truth`. The correct rejection (`Cor Rej`) and `Hit` are the right decisions. `Miss` and false alarm (`False Al`) are wrong decisions. The centre shows the cumulative tree performance in terms of mean of used cues (`mcu`), Percent of ignored cues (`pci`), sensitivity (`sens`), specificity (`spec`), accuracy (`acc`) and weighted Accuracy (`wacc`). The receiver operating characteristic (ROC) curve on the right-hand side compares the performance of all trees in the FFTrees object. The system also displays the performance of the fast frugal trees (`#`, green), CART (`C`, red), logistical regression (`L`, blue), random forest (`R`, violet) and the support vector machine (`S`, yellow).", fig.height=11, fig.width=11----
plot(res_mblrr.fft, decision.lables = c("Positive", "Negative"))

## ---- echo=FALSE, fig.cap="\\label{plot_peaks_ratio}Principle behind the `peaks_ratio` feature. The computation is based on a sequential linking of functions. The ``diffQ()`` function (\\texttt{MBmca}) determines numerically the first derivative of an amplification curve. This derivative is passed to the ``mcaPeaks()`` function (\\texttt{MBmca}). In the output all minima and all maxima are contained. The ranges are calculated from the minima and maxima. The Lagged Difference is determined from the ranges of the minima and maxima. Finally, the ratio of the differences (maximum/minimum) is calculated.", fig.height=8----

par(mfrow = c(2,1))

for(i in 1:2){
dat_smoothed <- chipPCR::smoother(RAS002[, 1], RAS002[, i+1])
res_diffQ <- suppressMessages(MBmca::diffQ(cbind(RAS002[-c(1,2), 1], dat_smoothed[-c(1,2)]), verbose = TRUE)$xy)
res_mcaPeaks <- MBmca::mcaPeaks(res_diffQ[, 1], res_diffQ[, 2])

range_p.max <- range(res_mcaPeaks$p.max[, 2])
diff_range_p.max <- diff(range_p.max)


range_p.min <- range(res_mcaPeaks$p.min[, 2])
diff_range_p.min <- diff(range_p.min)


peaks_ratio <- diff_range_p.max / diff_range_p.min


plot(res_diffQ, xlab = "Cycle", ylab = "dRFU/dCycle", type = "b", 
    ylim = c(diff_range_p.min, diff_range_p.max))

points(res_mcaPeaks$p.max, col = "red", pch = 19)
points(res_mcaPeaks$p.min, col = "cyan", pch = 19)

arrows(mean(res_mcaPeaks$p.max[, 1]), 
       range(res_mcaPeaks$p.max[, 2])[[1]], 
       mean(res_mcaPeaks$p.max[, 1]), 
       range(res_mcaPeaks$p.max[, 2])[[2]], col = "red", angle = 90, code = 3)

points(mean(res_mcaPeaks$p.max[, 1]), diff_range_p.max, col = "red", pch = 13, cex = 2)

arrows(mean(res_mcaPeaks$p.min[, 1]), 
       range(res_mcaPeaks$p.min[, 2])[[1]], 
       mean(res_mcaPeaks$p.min[, 1]), 
       range(res_mcaPeaks$p.min[, 2])[[2]], col = "cyan", angle = 90, code = 3)
       
points(mean(res_mcaPeaks$p.min[, 1]), diff_range_p.min, col = "cyan", pch = 13, cex = 2)

mtext(LETTERS[i], cex = 1.2, side = 3, adj = 0, font = 2)


legend("bottomleft", c("Local maxima", "Mean local maxima", 
"Local minima", "Mean local minima", paste("peaks_ratio", signif(diff_range_p.max/diff_range_p.min, 2))), 
                       col = c("red", "cyan", "red", "cyan", "black"), pch = c(19, 13, 19, 13, 1), 
                       bty = "n")
}

## ---- echo=FALSE, fig.cap="\\label{plot_cp_area}Analysis of area and changepoint features. Amplification curves from the datasets `stepone_std`, `RAS002`, `RAS003`, `lc96_bACTXY`, `C126EG595` and `dil4reps94` were analyzed with the ``encu()`` function. These datasets contain positive and negative amplification curves. Furthermore, the meta dataset contains amplification curves that exhibit a hook effect or non-sigmoid shapes, for instance. All amplification curves are manually classified. Altogether 626 positive and 317 negative amplification curves were included in the analysis. A) `polyarea`, is the area under the amplificatin curve determined by the Gauss polygon area formula. B) `peaks_ratio`, is the ratio of the local minima and the local maxima. C) `changepoint_e.agglo`, makes use of energy agglomerative clustering. Positive amplification curves have fewer change points than negative amplification curves. These two change point analyses generally separate positive and negative amplification curves. D) `changepoint_bcp`, analyses change points by a Bayesian approach. Positive amplification curves appear to contain more change points than negative amplification curves. Nevertheless, there is an overlap between the positive and negative amplification curves in both methods. This can lead to false-positive or false-negative classifications. E) `amptester_polygon`, is the cycle normalized order of a polygon.  F) `amptester_slope.ratio`, is the slop (linear model) of the raw fluorescence values at the approximate first derivate maximum, second derivative minimum and second derivative maximum.", fig.height=11, fig.width=11----
library(PCRedux)

data <- data_sample[data_sample$dataset %in% c("stepone_std", "RAS002", "RAS003", "lc96_bACTXY", "C126EG595", "dil4reps94"), ]

feature <- c(
    "polyarea", "peaks_ratio", "changepoint_e.agglo", "changepoint_bcp", 
  "amptester_polygon", "amptester_slope.ratio")

par(mfrow = c(2,3))

x <- data$decision

for(i in 1L:length(feature)) {
    y <- data[, colnames(data) == feature[i]]
    res <- stats::wilcox.test(y ~ x)
    h <- max(na.omit(y))
    l <- min(na.omit(y))
    h_text <- rep(h * 0.976, 2)
    
    par(bg=NA)
    stripchart(y ~ x, vertical = TRUE, ylab = feature[i],
               method = "jitter", pch = 19, cex = 1, 
               col = adjustcolor("darkgrey", alpha.f = 0.65), 
               ylim = c(l * 0.95, h * 1.05))
    
    boxplot(y ~ x, outline = FALSE, add = TRUE, boxwex = 0.35)

    legend("topleft", paste0("P = ", signif(res[["p.value"]])), 
           cex = 1, bty = "n")
    
    mtext(paste0(LETTERS[i], "   ", feature[i]), cex = 1.2, side = 3, 
          adj = 0, font = 2, col = ifelse(signif(res[["p.value"]]) < 0.05, 
                                          "black", "red"))
}

## ---- echo=FALSE, fig.cap="\\label{plot_cpa}Application of Bayesian change point analysis and energy agglomerative change point analysis methods to the `RAS002` dataset. An analysis of a negative and a positive amplification curve from the `RAS002` dataset was performed using the ``pcrfit_single()`` function. In this process, the amplification curves were analysed for change points using Bayesian change point analysis and energy agglomerative clustering. A) The negative amplification curve has a base signal of cica 2450 RFU and only a small signal increase to 2650 RFU. There is a clear indication of the signal variation (noise). B) The first negative derivative amplifies the noise so that some peaks are visible. C) The change point analysis shows changes in energy agglomerative clustering at several positions (green vertical line). The Bayesian change point analysis rarely exceeds a probability of 0.6 (grey vert line). D) The positive amplification curve has a lower base signal (~ 2450 RFU) and increases up to the 40th cycle (~3400 RFU). A sigmoid shape of the curve is clearly visible. E) The first negative derivation of the positive amplification curve shows a distinctive peak with a minimum at cycle 25. F) The change point analysis in energy agglomerative clustering shows changes (green vertical line) only at two positions. The Bayesian change point analysis shows a probability higher then 0.6 (grey horizontal line) at several positions.", fig.height=5----
options(warn = -1)

# Analyze a positiv and a negative amplification curve from the `RAS002` dataset 
# for change points using the `bcp` and `ecp` packages.
# 
library(bcp)
library(ecp)

# The MBmca package is used to calculate the approximate first derivative 
# of the amplification curve.
library(MBmca)

index <- which(grepl("B.Globin", colnames(RAS002)))

data <- RAS002[, c(1, index)]

amp_data <- data[, c(1,25, 2)]

# Smooth data with moving average for other data
# analysis steps.
dat_smoothed <- cbind(
  chipPCR::smoother(amp_data[, 1], amp_data[, 2]),
  chipPCR::smoother(amp_data[, 1], amp_data[, 3])
)

# Calculate the first derivative
dat_smoothed_deriv <- cbind(
    suppressMessages(MBmca::diffQ(cbind(amp_data[-c(1,2), 1], dat_smoothed[-c(1,2), 1]), verbose = TRUE)$xy)[, 2],
    suppressMessages(MBmca::diffQ(cbind(amp_data[-c(1,2), 1], dat_smoothed[-c(1,2), 2]), verbose = TRUE)$xy)[, 2]
)


# Bayesian analysis of change points
# Positive amplification curve
res_bcp_pos <- bcp(dat_smoothed_deriv[, 1])

y2_range <- range(res_bcp_pos$posterior.prob, na.rm = TRUE)

# Negative amplification curve
res_bcp_neg <- bcp(dat_smoothed_deriv[, 2])


# Change point analysis by energy agglomerative clustering

# Positive amplification curve
res_ecp_pos <- ecp::e.agglo(as.matrix(dat_smoothed_deriv[, 1]))$estimates
# Negative amplification curve
res_ecp_neg <- ecp::e.agglo(as.matrix(dat_smoothed_deriv[, 2]))$estimates

par(mfrow = c(2, 3))

plot(amp_data[, 1], amp_data[, 2], xlab = "Cycle", ylab = "RFU", type = "l", lwd = 2)
mtext("A    Negative", cex = 1.2, side = 3, adj = 0, font = 2)

plot(res_bcp_pos$data, xlab = "Cycle", ylab = "d(RFU) / d(cycle)", type = "l", lwd = 2)
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)

plot(res_bcp_pos$posterior.prob, xlab = "Cycle", ylab = "Probability", type = "b", lwd = 2, pch = 19, col = "red", ylim = c(0, 1))
res_ecp_pos <- ecp::e.agglo(as.matrix(dat_smoothed_deriv[, 1]))$estimates
abline(v = res_ecp_pos, col = "green", pch = 19)
abline(h = 0.6, col = "grey")
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2)
legend("topright", c("Baysian", "Agglomerative"), pch = c(19, 19), 
       col = c("red", "green"), bty="n")

plot(amp_data[, 1], amp_data[, 3], xlab = "Cycle", ylab = "RFU", type = "l", lwd = 2)
mtext("D    Positive", cex = 1.2, side = 3, adj = 0, font = 2)

plot(res_bcp_neg$data, xlab = "Cycle", ylab = "d(RFU) / d(cycle)", type = "l", lwd = 2)
mtext("E", cex = 1.2, side = 3, adj = 0, font = 2)

plot(res_bcp_neg$posterior.prob, xlab = "Cycle", ylab = "Probability", type = "b", lwd = 2, pch = 19, col = "red", ylim = c(0, 1))
res_ecp_neg <- ecp::e.agglo(as.matrix(dat_smoothed_deriv[, 2]))$estimates
abline(v = res_ecp_neg, col = "green", pch = 19)
abline(h = 0.6, col = "grey")
mtext("F", cex = 1.2, side = 3, adj = 0, font = 2)
legend("topleft", c("Baysian", "Agglomerative"), pch = c(19, 19), 
       col = c("red", "green"), bty="n")

## ------------------------------------------------------------------------
options(warn = -1)
suppressMessages(library(randomForest))
library(PCRedux)

data <- data_sample[data_sample$dataset %in% 
                c("batsch1",
                  "HCU32_aggR",
                  "lc96_bACTXY",
                  "RAS002", 
                  "RAS003", 
                  "stepone_std"), ]

n_positive <- sum(data[["decision"]] == "y")
n_negative <- sum(data[["decision"]] == "n")

dat <- data.frame(data[, c("amptester_shapiro", 
                           "amptester_lrt", 
                           "amptester_rgt", 
                           "amptester_tht", 
                           "amptester_slt",
                           "amptester_polygon", 
                           "amptester_slope.ratio")],
                  decision = as.numeric(factor(data$decision, 
                                         levels = c("n", "y"), 
                                         label = c(0, 1))) - 1)

# Select randomly observations from 70% of the data for training.
# n_train is the number of observations used for training.

n_train <- round(nrow(data) * 0.7)

# index_test is the index of observations to be selected for the training
index_test <- sample(1L:nrow(dat), size = n_train)

# index_test is the index of observations to be selected for the testing
index_training <- which(!(1L:nrow(dat) %in% index_test))

# train_data contains the data used for training

train_data <- dat[index_test, ]

# test_data contains the data used for training

test_data <- dat[index_training, ]


model_rf = randomForest(decision ~ ., data = train_data, ntree = 4000, 
                        importance = TRUE)


# Determine variable importance
res_importance <- importance(model_rf)

## ---- echo=TRUE, fig.cap="\\label{plot_random_forest}Random Forest.", fig.height=3----
par(mfrow = c(1,3))

plot(model_rf, main = "", las = 2)
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)

rownames(res_importance) <- substr(rownames(res_importance), 11, 22)


barplot(t(as.matrix(sort(res_importance[, 1]))), 
        ylab = "%IncMSE", main = "", las = 2,
        col = adjustcolor("grey", alpha.f = 0.5))
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)

barplot(t(as.matrix(sort(res_importance[, 2]))), 
        ylab = "IncNodePurity", main = "", las = 2,
        col = adjustcolor("grey", alpha.f = 0.5))
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2, las = 0)

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  # Copy and paste the code to an R console to evaluate it
#  
#  library(parallel)
#  library(doParallel)
#  
#  pcrfit_parallel <- function(data, n_cores=1) {
#      # Determine the number of available cores and register them
#      if (n_cores == "all") {
#          n_cores <- detectCores()
#      }
#  
#      registerDoParallel(n_cores)
#  
#      # Prepare the data for further processing
#      # Normalize RFU values to the alpha percentile (0.99)
#      cycles <- data.frame(cycles = data[, 1])
#      data_RFU <- data.frame(data[, -1])
#      data_RFU_colnames <- colnames(data_RFU)
#      data_RFU <- sapply(1L:ncol(data_RFU), function(i) {
#          data_RFU[, i] / quantile(data_RFU[, i], 0.99, na.rm = TRUE)
#      })
#      colnames(data_RFU) <- data_RFU_colnames
#  
#      # just to shut RCHeck for NSE we define ith_cycle
#      ith_cycle <- 1
#  
#      run_res <- foreach::foreach(
#          ith_cycle = 1L:ncol(data_RFU),
#                                  .packages = c(
#                                      "bcp", "changepoint", "chipPCR", "ecp", "MBmca",
#                                      "PCRedux", "pracma", "qpcR", "robustbase",
#                                      "zoo"
#                                  ),
#                                  .combine = rbind
#      ) %dopar% {
#          suppressMessages(pcrfit_single(data_RFU[, ith_cycle]))
#      }
#  
#  
#      res <- cbind(runs = colnames(data_RFU), run_res)
#  
#      rownames(res) <- NULL
#  
#      res
#  }
#  
#  # Calculate curve features of an amplification curve data. Note: Not all
#  # CPU cores are used. If need set "all" to use all available cores.
#  # In this example the testdat dataset from the qpcR package is used.
#  # The observations F1.1 and F1.2 are positive amplification curves. The observations
#  # F1.3 and F1.4 are negative.
#  options(warn = -1)
#  library(qpcR)
#  res_pcrfit_parallel <- pcrfit_parallel(testdat[, 1:5])
#  res_pcrfit_parallel

