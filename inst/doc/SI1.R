## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  # Select your local mirror
#  install.packages("PCRedux")

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  # The following command points to the help for download and install of packages
#  # from CRAN-like repositories or from local files.
#  ?install.packages()

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  # Install devtools, if you haven't already.
#  install.packages("devtools")
#  
#  library(devtools)
#  install_github("devSJR/PCRedux")

## ----eval=TRUE, echo=FALSE, results="asis"-------------------------------
# Load the RDML package for reading of the hookreg.rdml-file from the PCRedux
# package. The magrittr package is used for pipes.

library(RDML)
library(PCRedux)
library(magrittr)

# A comprehensive description of the RDML-file import can be found in Rödiger 
# et al. (2017) Bioinformatics

raw_data <- RDML$new(filename=system.file("hookreg.rdml", package="PCRedux"))

raw_data_tab <- raw_data$AsTable(name.pattern=paste(
    react$position,
    react$sample$id,
    # run id added to names
    sep="~"))
data <- as.data.frame(raw_data$GetFData(raw_data_tab, long.table=FALSE))

## ----ampcurves,fig.width=7.25, fig.height=9.25, echo=TRUE, message=FALSE, warning=FALSE, fig.cap="Amplification curves. A) Synthetic template, detected with Syto-13. B) Human \\textit{MLC-2v}, detected with a hydrolysis probe. C) \\textit{S27a} housekeeping gene, detected with SybrGreen I. D) Whole genome amplification, detected with EvaGreen. E) Human \\textit{BRCA1} gene, detected with a hydrolysis probe. F) Human \\textit{NRAS} gene, detected with a hydrolysis probe. G) Water control, detected with a hydrolysis probe. See Table \\ref{Table_human_rated} for details. RFU, relative fluorescence units. \\label{ampcurves}"----
par(mfrow=c(4,2))

# Plot all data of the hookreg.rdml-file according to their type.
# Synthetic template, detected with Syto-13
matplot(data[, 1], data[, 2:13], type="l", lty=1, lwd=2, ylab="RFU", xlab="Cycle")
mtext("A", cex = 1.8, side = 3, adj = 0, font = 2)

# Human MLC-2v, detected with a hydrolysis probe.
matplot(data[, 1], data[, 14:45], type="l", lty=1, lwd=2, ylab="RFU", xlab="Cycle")
mtext("B", cex = 1.8, side = 3, adj = 0, font = 2)

# S27a housekeeping gene, detected with SybrGreen I.
matplot(data[, 1], data[, 46:69], type="l", lty=1, lwd=2, ylab="RFU", xlab="Cycle")
mtext("C", cex = 1.8, side = 3, adj = 0, font = 2)

# Whole genome amplification, detected with EvaGreen.
matplot(data[, 1], data[, 70:71], type="l", lty=1, lwd=2, ylab="RFU", xlab="Cycle")
mtext("D", cex = 1.8, side = 3, adj = 0, font = 2)

# Human BRCA1 gene, detected with a hydrolysis probe.
matplot(data[, 1], data[, 72:87], type="l", lty=1, lwd=2, ylab="RFU", xlab="Cycle")
mtext("E", cex = 1.8, side = 3, adj = 0, font = 2)

# Human NRAS gene, detected with a hydrolysis probe.
matplot(data[, 1], data[, 88:95], type="l", lty=1, lwd=2, ylab="RFU", xlab="Cycle")
mtext("F", cex = 1.8, side = 3, adj = 0, font = 2)

# Water control, detected with a hydrolysis probe.
matplot(data[, 1], data[, 96:97], type="l", lty=1, lwd=2, ylab="RFU", xlab="Cycle")
mtext("G", cex = 1.8, side = 3, adj = 0, font = 2)

## ----eval=TRUE, echo=FALSE, results="asis"-------------------------------
library(readxl)
library(xtable)
options(xtable.comment=FALSE)

Table_human_rated <- read_xlsx(path=system.file("Table_human_rated.xlsx", 
                                                package="PCRedux"))

print(xtable(Table_human_rated, digits=0, 
             caption = "Overview of the used amplification curve data. The 
samples names, data source (origin of data either from an existing data set or  
prepared for this study), the detection chemistries (intercalator (Syto-13, 
SyberGreenI, EvaGreen), hydrolysis probes (TaqMan (Cy5/BHQ2) , TaqMan 
(HEX/BHQ1))) and calculations by tow humans.", 
             label='Table_human_rated'),
      size = "\\tiny",
      include.rownames = FALSE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.65'
)

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
# Load PCRedux package to obtain the data and make the hookreg() function
# available.
library(PCRedux)
# Load the magrittr to use the %>% pipe-operator
library(magrittr)

# `data` is a temporary data frame of the hook.rdml amplification curve data file.
# Apply the hookreg() function over the amplification curves and arrange the 
# results in the data frame `res_hookreg`.

res_hookreg <- sapply(2L:ncol(data), function(i) {
    hookreg(x=data[, 1], y=data[, i])
}) %>% t %>% data.frame(sample=colnames(data)[-1],.)

# Fetch the calculated parameters from the calculations with the hookreg() 
# function as a table `res_hookreg_table`.

res_hookreg_table <- data.frame(sample=as.character(res_hookreg[["sample"]]),
                                intercept=signif(res_hookreg[["intercept"]], 2),
                                slope=signif(res_hookreg[["slope"]], 1),
                                hook.start=signif(res_hookreg[["hook.start"]], 0),
                                hook.delta=signif(res_hookreg[["hook.delta"]], 0),
                                p.value=signif(res_hookreg[["p.value"]], 4),
                                CI.low=signif(res_hookreg[["CI.low"]], 2),
                                CI.up=signif(res_hookreg[["CI.up"]], 2),
                                hook.fit=res_hookreg[["hook.fit"]],
                                hook.CI=res_hookreg[["hook.CI"]],
                                hook=res_hookreg[["hook"]]
)

## ----results='asis', echo=TRUE, message=FALSE, warning=FALSE-------------
# Load the xtable to create a LaTeX table from the `res_hookreg_table`.
library(xtable)
options(xtable.comment=FALSE)
print(xtable(res_hookreg_table, 
             caption = "Results from the hookreg() function for the hookreg.rdml 
             data set.", 
             label='res_hookreg_table'),
      size = "\\tiny",
      include.rownames = FALSE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.65'
)

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
# Note that the PCRedux package and the magrittr package need to be loaded (see above).
# Load the qpcR package to prevent messages during the start.
suppressMessages(library(qpcR))

# `data` is a temporary data frame of the hook.rdml amplification curve data file.
# Apply the hookregNL() function over the amplification curves and arrange the 
# results in the data frame `res_hookregNL`.
# Not that `suppressMessages()` to prevent warning messages from the qpcR package.

res_hookregNL <- suppressMessages(sapply(2L:ncol(data), function(i) {
    hookregNL(x=data[, 1], y=data[, i])
}) %>% t %>% data.frame(sample=colnames(data)[-1],.))

res_hookregNL_table <- data.frame(sample=as.character(res_hookregNL[["sample"]]),
                                  slope=signif(as.numeric(res_hookregNL[["slope"]]), 1),
                                  CI.low=signif(as.numeric(res_hookregNL[["CI.low"]]), 2),
                                  CI.up=signif(as.numeric(res_hookregNL[["CI.up"]]), 2),
                                  hook.CI=unlist(res_hookregNL[["hook"]])
)

## ----results='asis', echo=TRUE, message=FALSE, warning=FALSE-------------
library(xtable)
options(xtable.comment=FALSE)

print(xtable(res_hookregNL_table, 
             caption = "Results from the hookregNL() function for the 
             hookreg.rdml data set.", 
             label='res_hookregNL_table'),
      size = "\\tiny",
      include.rownames = FALSE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.65'
)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
# Load the readxl package to obtain the classifications by the human.
# The classification data are stored in an EXCEL file, which is contained in the 
# PCRedux package.

library(readxl)

# The dplyr package was used for data manipulation.
library(dplyr)

# Read the data as EXCEL file.
Table_human_rated <- read_xlsx(path=system.file("Table_human_rated.xlsx", 
                                                package="PCRedux"))

# Aggregate the results from the statistical analysis and the human classifications.

res <- data.frame(sample=Table_human_rated[, "Sample"],
                  hr=Table_human_rated[, "Hook effect-like\nRater 2"],
                  hookreg=res_hookreg_table[, "hook"],
                  hookregNL=data.frame(unlist(res_hookregNL[, "hook"]))
)
colnames(res) <- c("Sample",
                   "Human rater",
                   "hookreg",
                   "hookregNL"
)

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
# A simple logic was applied to improve the classification result. In this case
# the assumption was, that an amplification curve has an hook effect or hook effect-like
# curvature, if either the hookreg() or hookregNL() function are positive.

meta_hookreg <- sapply(1:nrow(res), function(i){
    ifelse(res[i, "hookreg"] == 1 || res[i, "hookregNL"] == 1, 1, 0)
})

res_out <- data.frame(Sample=res[["Sample"]], res[["Human rater"]], 
                      res_hookreg[["hook"]], res_hookregNL_table[["hook.CI"]], 
                      meta_hookreg)

colnames(res_out) <- c("Sample",
                       "Human rater",
                       "hookreg",
                       "hookregNL",
                       "hookreg and hoohkreNL combined"
)

## ----results='asis', echo=TRUE, message=FALSE, warning=FALSE-------------
library(xtable)
options(xtable.comment=FALSE)

print(xtable(res_out, digits=0, 
             caption = "Aggregated decisions from the human classification and 
the results from the machine decision of the hookreg() and hookregNL() 
functions.", label='method_comparision'), ,
      caption.placement = "top",
      scalebox='0.65')

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
res_performeR <- rbind(
    hookreg=performeR(res_out[["hookreg"]], res_out[["Human rater"]]),
                       hookregNL=performeR(res_out[["hookregNL"]], res_out[["Human rater"]]),
                       combined_hookreg=performeR(res_out[["hookreg and hoohkreNL combined"]], 
                                                  res_out[["Human rater"]])
) %>% t %>% signif(4)

colnames(res_performeR) <- c("hookreg", "hookregNL", "hookreg and hookregNL")

## ----results='asis', echo=TRUE, message=FALSE, warning=FALSE-------------
library(xtable)
options(xtable.comment=FALSE)

print(xtable(res_performeR, digits=4, 
             caption = "Analysis of the performance of both algorithms. The 
performance of the individual test and the combination of the tests is shown. 
Note that the classification improved if the hookreg() and hookregNL() function 
were combined by a logical statement. The measure were determined with the 
\\textit{performeR()} function from the \\texttt{PCRedux} package. Sensitivity, 
TPR; Specificity, SPC; Precision, PPV; Negative predictive value, NPV; Fall-out, 
FPR; False negative rate, FNR; False discovery rate, FDR; Accuracy, ACC; F1 
score, F1; Matthews correlation coefficient, MCC, Cohen's kappa (binary 
classification), $\\kappa$", label='res_performeR'),
      size = "normalsize",
      include.rownames = TRUE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.75'
)

