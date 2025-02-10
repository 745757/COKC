# Load required libraries
library(GSVA)       # For gene set variation analysis (GSVA)
library(survival)   # For survival analysis
library(survminer)  # For survival curve visualization
library(survMisc)   # For optimal cutoff determination

# Step 1: Compute GSVA Score
# TPM is the TCGA-HNSC dataset in TPM format, containing survival information
list <- list()
list[["C3signature"]] <- gene  # Define the gene signature set
dat <- t(TPM)  # Transpose TPM for GSVA input

GSVA_C3_score <- gsva(expr = as.matrix(dat), 
                      gset.idx.list = list, 
                      mx.diff = FALSE, 
                      kcdf = "Gaussian", 
                      parallel.sz = 16) 

# Convert GSVA score matrix to dataframe for analysis
GSVA_C3_score <- as.data.frame(t(GSVA_C3_score))
colnames(GSVA_C3_score) <- "C3signature_score"

# Merge GSVA scores with clinical survival data
data <- merge(data, GSVA_C3_score, by = "row.names", all.x = TRUE)
colnames(data)[1] <- "SampleID"

# Step 2: Determine the Optimal Cutoff for C3signature Score
# Use maximally selected rank statistics to determine the best cutoff
cut <- surv_cutpoint(data, time = "OS", event = "EVENT", variables = "C3signature_score")
best_cutoff <- cut$cutpoint[["C3signature_score"]]  # Extract optimal cutoff value

# Assign high and low groups based on the optimal cutoff
data$group <- ifelse(data$C3signature_score > best_cutoff, "High", "Low")
data$group <- factor(data$group, levels = c("Low", "High"))

# Step 3: Perform Survival Analysis
# Conduct survival difference test
fitd <- survdiff(Surv(OS, EVENT) ~ group, data = data, na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# Fit survival curves
fit <- survfit(Surv(OS, EVENT) ~ group, data = data)

# Generate p-value label
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ", round(pValue, 3))))

# Step 4: Plot the Survival Curve
plot <- ggsurvplot(fit,
                   data = data,
                   pval = p.lab,
                   conf.int = TRUE,
                   risk.table = FALSE,
                   risk.table.col = "strata",
                   palette = c("#766DA7","#D79865"),
                   legend.labs = c("Low", "High"),
                   size = 1,
                   xlim = c(0, 120),
                   break.time.by = 20,
                   legend.title = "C3signature Score",
                   surv.median.line = "hv",
                   ylab = "Survival Probability (%)",
                   xlab = "Time (Months)",
                   ncensor.plot = FALSE,
                   ncensor.plot.height = 0.25,
                   risk.table.y.text = FALSE)

# Return the survival plot
return(plot)



















