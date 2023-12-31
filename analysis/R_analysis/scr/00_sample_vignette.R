library(RAPToR)
# vignette("RAPToR")

# RAPToR is a computational framework to estimate the real age of biological samples from gene expression. This is especially useful for fast-developing organisms – such as C. elegans worms, flies, or zebrafish – where many factors substantially impact developmental speed, and thus unintended developmental variation between samples could obscure or confound the effect of variables of interest.

# With RAPToR and the inferred age of your samples, you can

# precisely estimate the effect of perturbations on developmental timing,
# increase power in differential expression analyses,
# estimate differential expression due to uncontrolled development and,
# recover perturbation-specific effects on gene expression, even when the perturbation is completely confounded by development.


# Please cite our preprint (Bulteau and Francesconi (2021)) if you use RAPToR in your research:
# Bulteau R., Francesconi M. Real Age Prediction from the Transcriptome with RAPToR (2021) bioRxiv doi: 10.1101/2021.09.07.459270

# Why use RAPToR ?
# In gene expression data, unknown and unintended developmental variation among biological samples can obscure and confound the effect of variables of interest. As many factors can influence growth speed, synchronizing samples is particularly challenging, but failing to do so can strongly impact gene expression.

# Aware of the problem, studies with large scale developmental profiling generally re-order or rank the samples post-profiling with methods that combine dimension-reduction (e.g PCA, Diffusion Map) and a trajectory-finding method. Unfortunately, this is only feasible for experiments with hundreds of samples and/or time-series designs, which the overwhelming majority of expression profiling studies are not.

# The problem we are faced with is not limited to large scale experiments, so why should we only address it in such scenarios ?

# RAPToR provides a way to precisely determine the real age of single samples from their expression profile.

# How does it work ?
# The method works in a 2-step process.

# A reference gene expression time-series is interpolated to build a near-continuous, high-temporal-resolution reference (a number of which are included in associated data-packages, see below).
# A correlation profile against this reference is dressed for of each of your samples, and the timing of the correlation peak is the estimated age. Bootstrapping on genes then gives a confidence interval of the estimates.

# What type of data can be used ?
# RAPToR uses gene expression profiles to determine the age of your samples. This means that any method outputting information on gene expression on a large scale is appropriate : RNA-seq (preferably TPM), MicroArray…

# Note that the references provided in the data-packages are log(X+1) of expression values, so applying this transformation to your data is important when comparing expression changes with the reference (but is not required for staging).

# Warning :
#   Data must not be gene-centered, as this destroys the relationship between gene levels within a sample.


# General structure of RAPToR ---------------------------------------------
# The main package (RAPToR) holds all the necessary functions to stage samples and build references.

# Since we aimed to provide an easy way to predict the age of samples, we pre-built several references for commonly used organisms from available data in the literature. In R package standards, these are voluminous datasets so the references are stored in separate “data-packages.” For example, wormRef holds the C. elegans references.

# You can see the references available in a data-package with the list_refs() function (you must have the data-package installed for this). wormRef can be installed from this link.

list_refs("wormRef")

# The data-packages currently available are listed in the README of RAPToR’s github repo. You can also build your own reference data-packages following a few guidelines listed in the vignette on this topic.


# Usage example -----------------------------------------------------------
# In this part, we’ll show how to stage two C. elegans time-series datasets using references from wormRef (which you will need installed to use its references, see here for installation).

# A time-series of larval development in 4 different strains published by Aeschimann et al. (2017), hereafter called dsaeschimann2017.
# A high-resolution time-series of late larval development published by Hendriks et al. (2014), hereafter called dshendriks2014
# Both datasets are available on GEO (GSE80157, GSE52861).

# Loading the data
library(RAPToR)
# The code to create the dsaeschimann2017 and dshendriks2014 objects (downloading data) is lengthy and not the object of this vignette, you can find it at the end of the document.

# Here is what the data looks like

dsaeschimann2017$g[1:5,1:4]
#>                let.7.n2853._18hr let.7.n2853._20hr let.7.n2853._22hr let.7.n2853._24hr
#> WBGene00007063          7.501850         10.988212          10.45480          7.994587
#> WBGene00007064          8.023767          8.655388          14.21012          9.759401
#> WBGene00007065         15.919452         16.875057          15.23932         18.847718
#> WBGene00003525          1.416181         10.938876          13.42202          2.488798
#> WBGene00007067          1.765342          1.775650           2.77224          2.200257

head(dsaeschimann2017$p, n = 5)
#>                        title geo_accession           organism_ch1       strain
#> GSM2113587 let.7.n2853._18hr    GSM2113587 Caenorhabditis elegans let-7(n2853)
#> GSM2113588 let.7.n2853._20hr    GSM2113588 Caenorhabditis elegans let-7(n2853)
#> GSM2113589 let.7.n2853._22hr    GSM2113589 Caenorhabditis elegans let-7(n2853)
#> GSM2113590 let.7.n2853._24hr    GSM2113590 Caenorhabditis elegans let-7(n2853)
#> GSM2113591 let.7.n2853._26hr    GSM2113591 Caenorhabditis elegans let-7(n2853)
#>            time in development:ch1 age
#> GSM2113587                18 hours  18
#> GSM2113588                20 hours  20
#> GSM2113589                22 hours  22
#> GSM2113590                24 hours  24
#> GSM2113591                26 hours  26
dshendriks2014$g[1:5,1:5]
#>                contDevA_N2_21h contDevA_N2_22h contDevA_N2_23h contDevA_N2_24h contDevA_N2_25h
#> WBGene00007063       2.7298009        2.928176       2.3023107        2.960860        3.275716
#> WBGene00007064       5.2372118        6.161017       7.0817623        5.769718        6.910397
#> WBGene00007065      14.5928869       11.448793       9.7072795       11.616851       11.769170
#> WBGene00003525       0.2554172        1.182322       3.8705742       10.072059       16.148716
#> WBGene00007067       1.3860572        1.025319       0.6974292        1.031272        1.765380

head(dshendriks2014$p, n = 5)
#>                      title geo_accession time in development:ch1 age
#> GSM1277118 contDevA_N2_21h    GSM1277118                22 hours  22
#> GSM1277119 contDevA_N2_22h    GSM1277119                23 hours  23
#> GSM1277120 contDevA_N2_23h    GSM1277120                24 hours  24
#> GSM1277121 contDevA_N2_24h    GSM1277121                25 hours  25
#> GSM1277122 contDevA_N2_25h    GSM1277122                26 hours  26

# We first quantile-normalize and log the expression data.

dsaeschimann2017$g <- limma::normalizeBetweenArrays(dsaeschimann2017$g, method = "quantile")
dsaeschimann2017$g <- log1p(dsaeschimann2017$g) # log1p(x) = log(x + 1)

dshendriks2014$g <- limma::normalizeBetweenArrays(dshendriks2014$g, method = "quantile")
dshendriks2014$g <- log1p(dshendriks2014$g)

# You may also need to convert probe IDs or gene IDs of your data to match those of the reference series. For example, all the references included in the wormRef data-package use WormBase Gene IDs (e.g. WBGene00016153).

# We aggregate transcript-level data into gene-level expression in our references for a broader usage (since one can always compute gene-level expression from transcript expression, but not the other way around). RNASeq data should ideally be sum-aggregated at the count level, and other types of data – such as microarray, or RNASeq TPM (if counts are unavailable) – are mean-aggregated.

# To help with this conversion, gene ID reference tables are included in the data-packages with common gene IDs (built directly from biomaRt queries).

# Choosing a reference dataset
# RAPToR estimates the age of samples based on correlation with reference time-series. This means you have to select the proper reference to compare your samples with. You can determine which reference of a data-package is appropriate for your samples from the list with list_refs() like above, or by using the plot_refs() function.

# If we haven’t built references for your favorite organism yet, you can take a look at the Building your own references section in this vignette for a quick-start guide or the dedicated reference-building vignette for a more in-depth explanation.

# The chart below shows the references available in the wormRef data-package, along with landmark developmental stages.

plot_refs("wormRef")


# Loading the reference
# To increase the accuracy of age estimates beyond the temporal resolution of time-series profiling experiments, we interpolate on gene expression dynamics of the reference with respect to time. Pre-built references with optimal interpolation parameters can be loaded using the prepare_refdata() function.

# The Cel_larval reference would be appropriate for the two example datasets since the samples to stage are from mid-larval to early young adult. The n.inter parameter corresponds to the resolution of the interpolated reference. In the interest of lightening the computational load (each timepoint of the reference is compared to the samples) you can choose smaller values, but aim over 500 for optimal results.

r_larv <- prepare_refdata("Cel_larval", "wormRef", n.inter = 600)
# Note that age estimates will be given in the time unit and scale of the chosen reference (here, hours post-hatching at 20∘C).

# Age estimation
# All we need to do now is run the ae() (age estimation) function.

ae_dsaeschimann2017 <- ae(samp = dsaeschimann2017$g,                         # input gene expression matrix
                          refdata = r_larv$interpGE,            # reference gene expression matrix
                          ref.time_series = r_larv$time.series) # reference time-series
#> Bootstrap set size is 6239
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...
#> Warning in ae(samp = dsaeschimann2017$g, refdata = r_larv$interpGE, ref.time_series = r_larv$time.series): Some estimates come near the edges of the reference.
#> If possible, stage those on a different reference for confirmation.
#>                 nb.genes
#> refdata            18718
#> samp               19595
#> intersect.genes    18718
ae_dshendriks2014 <- ae(samp = dshendriks2014$g,                         # input gene expression matrix
                        refdata = r_larv$interpGE,            # reference gene expression matrix
                        ref.time_series = r_larv$time.series) # reference time-series
#> Bootstrap set size is 6239
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...
#> Warning in ae(samp = dshendriks2014$g, refdata = r_larv$interpGE, ref.time_series = r_larv$time.series): Some estimates come near the edges of the reference.
#> If possible, stage those on a different reference for confirmation.
#>                 nb.genes
#> refdata            18718
#> samp               19595
#> intersect.genes    18718
# Let’s look at the results.

plot(ae_dsaeschimann2017, groups = dsaeschimann2017$p$strain, show.boot_estimates = T)
plot(ae_dshendriks2014, show.boot_estimates = T)


# At 25∘C, C. elegans worms develop 1.5 times faster than worms at 20∘C. Since the worms we staged were grown at 25°C, but the reference – and thus, the estimated age – corresponds to 20∘C development, we can see this 1.5 factor by fitting a simple linear model between chronological and estimated age.

lm_dsaeschimann2017 <- lm(ae_dsaeschimann2017$age.estimates[,1] ~ dsaeschimann2017$p$age)
summary(lm_dsaeschimann2017)
#> 
#> Call:
#> lm(formula = ae_dsaeschimann2017$age.estimates[, 1] ~ dsaeschimann2017$p$age)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.06673 -0.64196  0.07043  0.59626  1.90407 
#> 
#> Coefficients:
#>                        Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)            -4.04996    0.66919  -6.052 3.65e-07 ***
#> dsaeschimann2017$p$age  1.50404    0.02352  63.950  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.9576 on 41 degrees of freedom
#> Multiple R-squared:  0.9901, Adjusted R-squared:  0.9898 
#> F-statistic:  4090 on 1 and 41 DF,  p-value: < 2.2e-16

lm_dshendriks2014 <- lm(ae_dshendriks2014$age.estimates[,1] ~ dshendriks2014$p$age)
summary(lm_dshendriks2014)
#> 
#> Call:
#> lm(formula = ae_dshendriks2014$age.estimates[, 1] ~ dshendriks2014$p$age)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -0.7566 -0.2409 -0.1657  0.1468  0.8074 
#> 
#> Coefficients:
#>                      Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)          -2.77413    0.74278  -3.735  0.00222 ** 
#> dshendriks2014$p$age  1.55958    0.02488  62.692  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.4587 on 14 degrees of freedom
#> Multiple R-squared:  0.9965, Adjusted R-squared:  0.9962 
#> F-statistic:  3930 on 1 and 14 DF,  p-value: < 2.2e-16

# Understanding the output ------------------------------------------------
# The output of ae() is an ae object with various elements, including age estimate and confidence intervals from booststrapping (age estimates on random gene subsets).


# General information can be accessed via the summary() function.

summary(ae_dshendriks2014)
#> 
#> Span of samples : 23.324
#> Range of samples :  [ 31.68 , 55.005 ]
#> -----------------------------------------------------------
#>                 age.estimate     lb     ub
#> contDevA_N2_21h       31.680 31.498 31.862
#> contDevA_N2_22h       33.150 32.968 33.332
#> contDevA_N2_23h       34.343 34.161 34.525
#> contDevA_N2_24h       35.904 35.859 35.950
#> contDevA_N2_25h       37.557 37.511 37.603
#> contDevA_N2_26h       39.210 39.096 39.324
#> contDevA_N2_27h       41.690 41.644 41.736
#> contDevA_N2_28h       43.159 43.045 43.273
#> contDevA_N2_29h       44.169 44.123 44.215
#> contDevA_N2_30h       45.363 45.249 45.477
#> contDevA_N2_31h       46.924 46.742 47.106
#> contDevA_N2_32h       48.485 48.439 48.531
#> contDevA_N2_33h       49.495 49.313 49.677
#> contDevA_N2_34h       51.423 51.241 51.605
#> contDevA_N2_35h       54.178 53.996 54.360
#> contDevA_N2_36h       55.005 54.959 55.050
#> -----------------------------------------------------------

# Age estimates and their confidence intervals are accessible through $age.estimates.
head(ae_dshendriks2014$age.estimates)
#>                 age.estimate       lb       ub cor.score
#> contDevA_N2_21h     31.68042 31.49836 31.86248 0.9312136
#> contDevA_N2_22h     33.14966 32.96760 33.33171 0.9310279
#> contDevA_N2_23h     34.34341 34.16135 34.52547 0.9346601
#> contDevA_N2_24h     35.90448 35.85856 35.95039 0.9328182
#> contDevA_N2_25h     37.55737 37.51145 37.60328 0.9330641
#> contDevA_N2_26h     39.21026 39.09627 39.32424 0.9359745

# The table holds the following :
#age.estimate, the global estimate for the sample (whole gene set).
# lb, ub, the lower and upper bounds of the bootstrapped age estimates’ confidence interval (Median Absolute Deviation).
# cor.score, the correlation score between the sample and reference at the age estimate.

# Plotting ----------------------------------------------------------------
# Estimates and confidence intervals can be displayed in the form of a ‘dotchart’ with the default plot() function (as done above).

# The ae object also holds record of correlation scores of each samples with the full reference span (for bootstrap estimates as well). These correlation profiles can be plotted with plot_cor.ae()

par(mfrow=c(2,2))
plot_cor.ae(ae_dshendriks2014, subset = c(1,4,9,11))
# The red bars correspond to the estimate confidence interval, and the 95% interval of bootstrap correlation with the reference is shown as black dotted lines. The sample age estimate is displayed below the interval.

# Building your own references --------------------------------------------
# This section is just a broad overview, there is a vignette entirely dedicated to reference-building.

# Interpolating on time-series expression data is the key to get the high-temporal-resolution references we need for RAPToR. When using our pre-built references this is done internally by calling the prepare_refdata() function, but what’s going on behind the scenes uses the functions described below.

# To build your reference, you will require a time-series of gene expression data for your favorite organism. The dsaeschimann2017 example dataset we loaded earlier will be used to illustrate the process.

# The gene expression interpolation model interface
# Gene expression interpolation models (GEIMs), are built with the ge_im() function. This function takes as input 3 key arguments :
  
# X : your time-series gene expression matrix (genes as rows, samples as columns)
# p : a dataframe of phenotypic data, samples as rows in the same order as X columns. This should include the age/time variable and any other covariate(s) you want to include in the model (e.g batch, strain)
# formula : the model formula. This should be a standard R formula using terms found in p, which may include elements (such as splines) from chosen model type (see below). It must start with X ~.
# For example, using the dsaeschimann2017 dataset we could build the following model.

m_dsaeschimann2017 <- ge_im(X = dsaeschimann2017$g, p = dsaeschimann2017$p, 
                            formula = "X ~ s(age, bs = 'ts') + strain", nc = 32)

# About the models --------------------------------------------------------
# In order to model such a large number of output variable (genes), our strategy is to project the data in a dimensionally-reduced space and interpolate there before re-projecting the data back to genes. We propose to do this on Principal Components or Independant Components ( Independant Component Analysis ).

# Both PCA and ICA perform the same type of linear transformation on the data (they just optimize different criteria). We get the following :
# X(m×n)=G(m×c)ST(n×c)
# with X, the matrix of m genes by n samples, G the gene loadings (m genes by c components) and ST the sample scores (n samples by c components). When performing PCA (or ICA) on gene expression data, S is what’s usually plotted (e.g. PC1 vs. PC2) to see how samples are grouped in the component space.

# Alter, Brown, and Botstein (2000) demonstrated that singular value decomposition of gene expression data can be taken as “eigengenes,” giving a global picture of the expression landscape and dynamics with a few components. We essentially use the same property for a GEIM. We fit a model on the columns of ST (eigengenes), predict in the component space, and reconstruct the gene expression data by a matrix product with the gene loadings.

# We’ve implemented 2 model types : Generalized Additive Models (GAMs, the default) and Generalized Linear Models (GLMs). GAMs rely on the gam() function of the mgcv package, and GLMs on the glm() function of the stats core package. The specified model formula can make use of all the tools one can use with gam() or glm(), most notably the variety of polynomial or smoothing splines implemented through the s() function for GAMs.

# If you are not familiar with the mgcv package or GLMs in R, we recommend you look at their function documentations (especially s() for GAMs) to understand what can be included in the formulas to handle non-linear dynamics.

# Note that a single model formula is specified and applied to all the components, but the models are fitted independently on the components.

# Finding the appropriate model and parameters
# Model type
# The default GEIM is fitted using GAMs on PCA components, which is a robust choice when applying a smoothing spline to the data. PCA and ICA interpolation yield near-identical results in most scenarios.
# 
# Parameter estimation
# The number of components to use for the interpolation is by default set to the number of samples. However, we recommend to set a cutoff on explained variance of PCA components to select it. For example, on the dsaeschimann2017 dataset, we set the threshold at 99% :
#   
#   pca_dsaeschimann2017 <- stats::prcomp(t(dsaeschimann2017$g), center = TRUE, scale = FALSE, rank = 25)
# nc <- sum(summary(pca_dsaeschimann2017)$importance[3,] < .99) + 1
# nc
# #> [1] 32
# Note that this threshold must be set with respect to the noise in the data. For example, in very noisy data, would you consider that 99% of the variance in the dataset corresponds to meaningful dynamics ? One can also keep only components that have ‘intelligible dynamics’ with respect to time, defined as those where a model fit explains >0.5 of the deviance.
# 
# 
# 
# Choosing from different splines (and/or parameters) can be done with Cross-Validation (CV) through the use of the ge_imCV() function. The function inputs the X, p and a formula_list to test. Other parameters on the CV itself can also be given (e.g. training set size).
# 
# Below is an example of usage to choose among available spline types for the dsaeschimann2017 GEIM.
# 
# smooth_methods <- c("tp", "ts", "cr", "ps")
# flist <- as.list(paste0("X ~ s(age, bs = \'", smooth_methods, "\') + strain")) #, k=", k, ", fx=TRUE
# flist
# #> [[1]]
# #> [1] "X ~ s(age, bs = 'tp') + strain"
# #> 
# #> [[2]]
# #> [1] "X ~ s(age, bs = 'ts') + strain"
# #> 
# #> [[3]]
# #> [1] "X ~ s(age, bs = 'cr') + strain"
# #> 
# #> [[4]]
# #> [1] "X ~ s(age, bs = 'ps') + strain"
# 
# cv_dsaeschimann2017 <- ge_imCV(X = dsaeschimann2017$g, p = dsaeschimann2017$p, formula_list = flist,
#                                cv.n = 20, nc = nc, nb.cores = 3)
# #> CV on 4 models. cv.n = 20 | cv.s = 0.8
# #> 
# #> ...Building training sets
# #> ...Setting up cluster
# #> ...Running CV
# #> ...Cleanup and formatting
# plot(cv_dsaeschimann2017, names = paste0("bs = ", k), outline = F,
#      swarmargs = list(cex = .8))
# 
# 
# ge_imCV() computes various indices of model performance : the average Correlation Coefficient (aCC), the average Relative Error (aRE), Mean Squared Error (MSE) and average Root MSE (aRMSE). These indices all compare model predictions and the true data. The ge_imCV() function computes them both on the validation set (CV Error) and on the training set (Model PerFormance).
# 
# From the plots above, we can see the different splines perform similarly in this scenario, all could work. We chose ts (a thin-plate regression spline), as it appears to minimize CV error without much impact model performance.
# 
# Predicting from the model
# The geim object returned by ge_im() has its predict() method which can be used like for any R model. You can also specify if you want the predictions in component space or as the full gene expression matrix.
# 
# On our dsaeschimann2017 example :
#   
#   # setup new data
#   n.inter <- 100
# ndat <- data.frame(age = seq(min(dsaeschimann2017$p$age), max(dsaeschimann2017$p$age),  l = n.inter), 
#                    strain = rep("N2", n.inter))
# 
# # predict
# pred_dsaeschimann2017_comp <- predict(m_dsaeschimann2017, ndat, as.c = TRUE) # in component space
# pred_dsaeschimann2017_ge <- predict(m_dsaeschimann2017, ndat)
# Checking/Validating the interpolation
# After building our model, we can look at the interpolation results by:
#   
#   Checking the model predictions against components (plots)
# Staging the samples on their own interpolated data, or better (if possible) stage another independent time-series on your reference for external validation.
# We can do both with our example, using the dshendriks2014 dataset for external validation.
# 
# 
# 
# # make a 'reference object' 
# r_dsaeschimann2017 <- list(interpGE = pred_dsaeschimann2017_ge, time.series = ndat$age)
# 
# ae_test_dsaeschimann2017 <- ae(dsaeschimann2017$g, r_dsaeschimann2017$interpGE, r_dsaeschimann2017$time.series)
# ae_test_dshendriks2014 <- ae(dshendriks2014$g, r_dsaeschimann2017$interpGE, r_dsaeschimann2017$time.series)
# 
# 
# 
# 
# Using a prior
# When few genes are available for staging (e.g. for tissue-specific staging), it may be appropriate to use a prior to help age estimation. In the ae() function, priors work by giving the parameters for gaussian distributions of time (for each sample). The correlation peaks are then ranked according to the prior’s density. The correlation profile is unaffected by the prior, only the choice of the correlation peak is.
# 
# This implies that with a prior which is completely off, the estimate may also be wrong ; use with care.
# 
# The priors are given in the reference series’ time scale, so beware of growth speed difference with temperature or different time origins (fertilization, egg-laying, hatching…). For example, the dshendriks2014 C. elegans data we used previously is grown at 25∘C, in contrast with the references’ time for 20∘C development.
# 
# Due to these possible differences and the bias introduced by the prior, we recommend to carefully plan its use. Performing a first run without priors will give a general idea of the difference between the chronological and developmental age of your samples.
# 
# Once the priors are determined, you will need to set the standard deviation of the gaussian centered on the sample with the prior.params argument. This parameter will also indirectly change the weight of the prior over the correlation score for estimate selection.
# 
# On our dshendriks2014 example, we can use ajusted known chronological ages for 20∘C.
# 
# priors <- dshendriks2014$p$age * 1.6 - 5 # rough approximation based on our previous lm
# 
# ae_dshendriks2014_prior <- ae(samp = dshendriks2014$g,
#                               refdata = r_larv$interpGE,
#                               ref.time_series = r_larv$time.series,
#                               prior = priors,
#                               prior.params = 10)
# #>                 nb.genes
# #> refdata            18718
# #> samp               19595
# #> intersect.genes    18718
# #> Bootstrap set size is 6239
# #> Performing age estimation...
# #> Bootstrapping...
# #>  Building gene subsets...
# #>  Computing correlations...
# #>  Performing age estimation...
# #> Computing summary statistics...
# #> Warning in ae(samp = dshendriks2014$g, refdata = r_larv$interpGE, ref.time_series = r_larv$time.series, : Some estimates come near the edges of the reference.
# #> If possible, stage those on a different reference for confirmation.
# plot(ae_dshendriks2014_prior, main="Age estimates with priors on dshendriks2014", show.boot_estimates = T,
#      show.prior = T, col.p = 'red', l.pos = 'bottomright')
# 
# all(ae_dshendriks2014_prior$age.estimates[,1]==ae_dshendriks2014$age.estimates[,1])
# #> [1] TRUE
# 
# 
# As you can see, here the estimates are identical to those without priors.
# 
# 
# Loading dsaeschimann2017 and dshendriks2014
# Note : set the data_folder variable to an existing path on your system where you want to store the objects.
# 
# data_folder <- "../inst/extdata/"
# 
# requireNamespace("wormRef", quietly = T)
# requireNamespace("utils", quietly = T)
# requireNamespace("GEOquery", quietly = T) # May need to be installed with bioconductor
# requireNamespace("Biobase", quietly = T) # same
# Utility functions
# 
# raw2tpm <- function(rawcounts, genelengths){
#   if(nrow(rawcounts) != length(genelengths))
#     stop("genelengths must match nrow(rawcounts).")
#   x <- rawcounts/genelengths
#   return(t( t(x) * 1e6 / colSums(x) ))
# }
# 
# fpkm2tpm <- function(fpkm){
#   return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
# }
# dsaeschimann2017
# 
# geo_dsaeschimann2017 <- "GSE80157"
# 
# g_url_dsaeschimann2017 <- GEOquery::getGEOSuppFiles(geo_dsaeschimann2017, makeDirectory = FALSE, fetch_files = FALSE)
# g_file_dsaeschimann2017 <- paste0(data_folder, "dsaeschimann2017.txt.gz")
# utils::download.file(url = as.character(g_url_dsaeschimann2017$url[2]), destfile = g_file_dsaeschimann2017)
# 
# X_dsaeschimann2017 <- read.table(gzfile(g_file_dsaeschimann2017), h=T, sep = '\t', stringsAsFactors = F, row.names = 1)
# 
# # convert to tpm & wb_id
# X_dsaeschimann2017 <- X_dsaeschimann2017[rownames(X_dsaeschimann2017)%in%wormRef::Cel_genes$wb_id,]
# X_dsaeschimann2017 <- raw2tpm(rawcounts = X_dsaeschimann2017, 
#                               genelengths = wormRef::Cel_genes$transcript_length[match(rownames(X_dsaeschimann2017),
#                                                                                        wormRef::Cel_genes$wb_id)])
# 
# # pheno data
# P_dsaeschimann2017 <- Biobase::pData(GEOquery::getGEO(geo_dsaeschimann2017, getGPL = F)[[1]])
# P_dsaeschimann2017[,10:34] <- NULL
# P_dsaeschimann2017[, 3:8] <- NULL
# 
# colnames(P_dsaeschimann2017)[4] <- "strain"
# P_dsaeschimann2017$strain <- factor(P_dsaeschimann2017$strain)
# P_dsaeschimann2017$title <- make.names(P_dsaeschimann2017$title)
# 
# colnames(X_dsaeschimann2017) <- gsub('RNASeq_riboM_', '', colnames(X_dsaeschimann2017), fixed = T)
# P_dsaeschimann2017$title <- gsub('RNASeq_riboM_', '', P_dsaeschimann2017$title, fixed = T)
# 
# # get age 
# P_dsaeschimann2017$age <- as.numeric(sub('(\\d+)\\shours', '\\1', P_dsaeschimann2017$`time in development:ch1`))
# 
# 
# X_dsaeschimann2017 <- X_dsaeschimann2017[, P_dsaeschimann2017$title]
# 
# dsaeschimann2017 <- list(g = X_dsaeschimann2017, p = P_dsaeschimann2017)
# save(dsaeschimann2017, file = paste0(data_folder, "dsaeschimann2017.RData"), compress = "xz")
# 
# # cleanup
# file.remove(g_file_dsaeschimann2017)
# rm(geo_dsaeschimann2017, g_url_dsaeschimann2017, g_file_dsaeschimann2017, X_dsaeschimann2017, P_dsaeschimann2017)
# dshendriks2014
# 
# geo_dshendriks2014 <- "GSE52861"
# 
# g_url_dshendriks2014 <- GEOquery::getGEOSuppFiles(geo_dshendriks2014, makeDirectory = FALSE, fetch_files = FALSE)
# g_file_dshendriks2014 <- paste0(data_folder, "dshendriks2014.txt.gz")
# utils::download.file(url = as.character(g_url_dshendriks2014$url[2]), destfile = g_file_dshendriks2014)
# 
# X_dshendriks2014 <- read.table(gzfile(g_file_dshendriks2014), h=T, sep = '\t', stringsAsFactors = F, row.names = 1)
# 
# # convert to tpm & wb_id
# X_dshendriks2014 <- X_dshendriks2014[rownames(X_dshendriks2014)%in%wormRef::Cel_genes$wb_id,]
# X_dshendriks2014 <- raw2tpm(rawcounts = X_dshendriks2014, 
#                             genelengths = wormRef::Cel_genes$transcript_length[match(rownames(X_dshendriks2014),
#                                                                                      wormRef::Cel_genes$wb_id)])
# 
# 
# # pheno data
# P_dshendriks2014 <- Biobase::pData(GEOquery::getGEO(geo_dshendriks2014, getGPL = F)[[1]])
# 
# # filter relevant fields/samples
# P_dshendriks2014 <- P_dshendriks2014[(P_dshendriks2014$`strain:ch1` == 'N2') & (P_dshendriks2014$`growth protocol:ch1` == 'Continuous'), ]
# P_dshendriks2014 <- P_dshendriks2014[, c("title", "geo_accession", "time in development:ch1")]
# 
# # get age 
# P_dshendriks2014$age <- as.numeric(sub('(\\d+)\\shours', '\\1', P_dshendriks2014$`time in development:ch1`))
# 
# 
# # formatting
# P_dshendriks2014$title <- gsub('RNASeq_polyA_', '', 
#                                gsub('hr', 'h', 
#                                     gsub('-', '.', fixed = T, as.character(P_dshendriks2014$title))))
# colnames(X_dshendriks2014) <- gsub('RNASeq_polyA_','', colnames(X_dshendriks2014))
# X_dshendriks2014 <- X_dshendriks2014[, P_dshendriks2014$title]
# 
# dshendriks2014 <- list(g = X_dshendriks2014, p = P_dshendriks2014)
# save(dshendriks2014, file = paste0(data_folder, "dshendriks2014.RData"), compress = "xz")
# 
# # cleanup
# file.remove(g_file_dshendriks2014)
# rm(geo_dshendriks2014, g_url_dshendriks2014, g_file_dshendriks2014, X_dshendriks2014, P_dshendriks2014)
# References
# Aeschimann, Florian, Pooja Kumari, Hrishikesh Bartake, Dimos Gaidatzis, Lan Xu, Rafal Ciosk, and Helge Großhans. 2017. “Lin41 Post-Transcriptionally Silences mRNAs by Two Distinct and Position-Dependent Mechanisms.” Molecular Cell 65 (3): 476–89.
# Alter, Orly, Patrick O Brown, and David Botstein. 2000. “Singular Value Decomposition for Genome-Wide Expression Data Processing and Modeling.” Proceedings of the National Academy of Sciences 97 (18): 10101–6.
# Bulteau, Romain, and Mirko Francesconi. 2021. “Real Age Prediction from the Transcriptome with RAPToR.” bioRxiv. https://doi.org/10.1101/2021.09.07.459270.
# Hendriks, Gert-Jan, Dimos Gaidatzis, Florian Aeschimann, and Helge Großhans. 2014. “Extensive Oscillatory Gene Expression During c. Elegans Larval Development.” Molecular Cell 53 (3): 380–92.