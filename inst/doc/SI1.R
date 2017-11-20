## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
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

# Plot all data of the hookreg.rdml-file
matplot(data[, 1], data[, -1], type="l", lty=1, lwd=2, ylab="RFU", xlab="Cycle")

## ----results='asis', echo=TRUE, message=FALSE, warning=FALSE-------------
library(xtable)
options(xtable.comment=FALSE)
print(xtable(t(rbind(head(data, 10), tail(data, 10)))), scalebox='0.55', floating=FALSE)

## ----eval=TRUE,echo=FALSE,results="asis"---------------------------------
library(readxl)
library(xtable)
options(xtable.comment=FALSE)

Table_human_rated <- read_xlsx(path=system.file("Table_human_rated.xlsx", package="PCRedux"))

print(xtable(Table_human_rated, digits=0),
      size = "normalsize",
      include.rownames = FALSE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.55', floating=FALSE
)

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
library(PCRedux)
library(magrittr)
suppressMessages(library(qpcR))

res_hookreg <- sapply(2L:ncol(data), function(i) {
  hookreg(x=data[, 1], y=data[, i])
}) %>% t %>% data.frame(sample=colnames(data)[-1],.)


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
library(xtable)
options(xtable.comment=FALSE)
print(xtable(res_hookreg_table),
      size = "normalsize",
      include.rownames = FALSE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.55', floating=FALSE
)

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
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

print(xtable(res_hookregNL_table),
      size = "normalsize",
      include.rownames = FALSE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.55', floating=FALSE
)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
library(readxl)
library(dplyr)

Table_human_rated <- read_xlsx(path=system.file("Table_human_rated.xlsx", package="PCRedux"))

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

## ----results='asis', echo=TRUE, message=FALSE, warning=FALSE-------------
library(xtable)
options(xtable.comment=FALSE)

print(xtable(res, digits=0),
      size = "normalsize",
      include.rownames = FALSE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.55', floating=FALSE
)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
meta_hookreg <- sapply(1:nrow(res), function(i){
  ifelse(res[i, "hookreg"] == 1 || res[i, "hookregNL"] == 1, 1, 0)
})

res_out <- data.frame(Sample=res[["Sample"]], res[["Human rater"]], res_hookreg[["hook"]], res_hookregNL_table[["hook.CI"]], meta_hookreg)

colnames(res_out) <- c("Sample",
                       "Human rater",
                       "hookreg",
                       "hookregNL",
                       "hookreg and hoohkreNL combined"
)

## ----results='asis', echo=TRUE, message=FALSE, warning=FALSE-------------
library(xtable)
options(xtable.comment=FALSE)

print(xtable(res_out, digits=0), scalebox='0.55', floating=FALSE)

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

print(xtable(res_performeR, digits=0),
      size = "normalsize",
      include.rownames = TRUE,
      include.colnames = TRUE,
      caption.placement = "top",
      comment=FALSE,
      table.placement = "!ht", scalebox='0.75', floating=FALSE
)

