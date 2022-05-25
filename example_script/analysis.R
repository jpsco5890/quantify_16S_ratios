# Load required libraries
library(CASEU) # Contains fitting functions
library(sangerseqR) # Reads '.ab1' files

# Read in the '.ab1' files mixtures are 1A01 to 3B05 (##_##)
pure_1A01_ab1 <- read.abif("sequences/1A01_100__3B05_000-27F-JR.ab1")
pure_3B05_ab1 <- read.abif("sequences/1A01_000__3B05_100-27F-JR.ab1")
mix_75_25_ab1 <- read.abif("sequences/1A01_075__3B05_025-27F-JR.ab1")
mix_50_50_ab1 <- read.abif("sequences/1A01_050__3B05_050-27F-JR.ab1")
mix_25_75_ab1 <- read.abif("sequences/1A01_025__3B05_075-27F-JR.ab1")

# Extract fluorescence signal for the four base-pair channels
# Construct a matrix from these vectors, and normalize to mean of 1
pure_1A01_nfm <- extractElectropherogram(pure_1A01_ab1,
                                         normLim = c(1500,9000))
pure_3B05_nfm <- extractElectropherogram(pure_3B05_ab1,
                                         normLim = c(1500,9000))
mix_75_25_nfm <- extractElectropherogram(mix_75_25_ab1,
                                         normLim = c(1500,9000))
mix_50_50_nfm <- extractElectropherogram(mix_50_50_ab1,
                                         normLim = c(1500,9000))
mix_25_75_nfm <- extractElectropherogram(mix_25_75_ab1,
                                         normLim = c(1500,9000))

# Plot the data to make sure it looks like an electropherogram
matplot(pure_1A01_nfm, 
        type = 'l', 
        lty = 1,
        xlab = "Time [indeces]")

# Zoom in to get a better look
matplot(pure_1A01_nfm, 
        type = 'l', 
        lty = 1,
        xlim = c(1000,2000),
        xlab = "Time [indeces]")

# Fit the mixtures to the pure electropherograms
pure_1A01_res <- fitSangerMixture(mixture = pure_1A01_nfm,
                                  components = list(s1A01 = pure_1A01_nfm,
                                                    s3B05 = pure_3B05_nfm),
                                  verbose = TRUE)
pure_3B05_res <- fitSangerMixture(mixture = pure_3B05_nfm,
                                  components = list(s1A01 = pure_1A01_nfm,
                                                    s3B05 = pure_3B05_nfm),
                                  verbose = TRUE)
mix_75_25_res <- fitSangerMixture(mixture = mix_75_25_nfm,
                                  components = list(s1A01 = pure_1A01_nfm,
                                                    s3B05 = pure_3B05_nfm),
                                  verbose = TRUE)
mix_50_50_res <- fitSangerMixture(mixture = mix_50_50_nfm,
                                  components = list(s1A01 = pure_1A01_nfm,
                                                    s3B05 = pure_3B05_nfm),
                                  verbose = TRUE)
mix_25_75_res <- fitSangerMixture(mixture = mix_25_75_nfm,
                                  components = list(s1A01 = pure_1A01_nfm,
                                                    s3B05 = pure_3B05_nfm),
                                  verbose = TRUE)

# Print the results to identify ratios and R-squared for the fit
print(pure_1A01_res)
print(pure_3B05_res)
print(mix_75_25_res)
print(mix_50_50_res)
print(mix_25_75_res)

# Manually inspect the fit
inspect_fit <- function(result = NULL, filename = NA){
  if(!is.na(filename)){
    png(filename = paste(filename, ".png", sep = ""),
        width = 900,
        height = 800)
  }
  ran_start <- sample(x = 2000:6000, size = 1)
  op <- par(mfrow=c(2,1), mar=c(4,4,1,1))
  plot(result)
  plot(result, xlim=c(ran_start,ran_start+500))
  par(op)
  if(!is.na(filename)){
    dev.off()
  }
}

inspect_fit(result = pure_1A01_res, filename = "pure_1A01_fit")
inspect_fit(result = pure_3B05_res, filename = "pure_3B05_fit")
inspect_fit(result = mix_75_25_res, filename = "mix_75_25_fit")
inspect_fit(result = mix_50_50_res, filename = "mix_50_50_fit")
inspect_fit(result = mix_25_75_res, filename = "mix_25_75_fit")

# Assess quality of fit with CASEU diagnostics
sangerFitDiagnosticPlot(pure_1A01_res)
sangerFitDiagnosticPlot(pure_3B05_res)
sangerFitDiagnosticPlot(mix_75_25_res)
sangerFitDiagnosticPlot(mix_50_50_res)
sangerFitDiagnosticPlot(mix_25_75_res)
