
# INDExA (Interval Null Differential EXpression/Abundance)

<!-- badges: start -->
<!-- badges: end -->

INDExA is a differential abundance/expression tool for the analysis
of DNA sequence count data (or potentially other sequence count data).
It is useful for analyzing RNA-seq or 16S rRNA-seq. It allows researchers
to define interval assumptions to account for uncertainty in normalizations
and/or the unmeasured scale. These interval assumptions are integrated into
the null hypothesis testing framework itself through interval null hypothesis
tests.

This is the alpha version of the INDExA package used for the analyses in the paper
"Replacing Normalizations with Interval Assumptions Improves the Rigor and Robustness
of Differential Expression and Differential Abudnance Analysis."

Currently this package supports simple analyses comparing two binary conditions with a
single interval assumption. We plan to extend this package to include linear models,
potentially including continous covariates and multiple interval support.

## Installation

You can install the development version of INDExA from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("Silverman-Lab/INDExA")
```
## Example

The following is a basic example analysis, the full analysis can be found in
the vignettes.

Here we analyze simulated data from a mock experiment measuring changes in
50 microbes' abundance in the colons of 100 individuals. The individuals were
divided into two groups: 50 *untreated* individuals representing control group
and 50 *treated* individuals who ingested a selective microbial growth factor.
30 of the microbes were simulated with positive LFC between 0.25 and 2.5, and
the remaining 20 microbes were simulated with LFC of 0.

Consider, based on prior knowledge, we believe the selective microbial growth
factor could only increase the total number of microbes in the colons of
**treated** vs. **untreated** individuals by at most 2.5 fold (and at least 1
fold). If we represent the $log_2$ fold change in total colonic microbial
abundance as the variable $\theta^\perp$, we can represent this (correct)
interval assumption as:

$$
\begin{align*}
  \theta^\perp &\in [\log_2(1), \log_2(2.5)] \\
  \theta^\perp &\in [0, 1.322].
\end{align*}
$$

(**NOTE:** In this sim the true $\theta^\perp=1.195$)

Using this interval assumption, we can perform an analysis using `INDExA`
function:

```R
library(INDExA)
data(data_sim_1)
# 0 represents untreated samples and 1 represents treated samples
conds <- c(rep(0, 50), rep(1, 50))
# Our interval assumption about the change in total microbial load
epsilon.interval <- c(0, 1.322)
# Build INDExA Object from MC Sampling
indexa.obj <- indexa.mc(data_sim_1$seq_counts, mc.samples=500)
# Run Interval Null Hypothesis Testing
indexa.sim.results <- indexa.test(indexa.obj, conds,
                                  epsilon.interval=epsilon.interval)
```
## INDExA Output Explained

The matrix returned from `indexa.test` has the following columns:

1. **median.effect**: The median LFC given the choice of normalization. NOTE: this does not account for the interval assumption
2. **median.effect.l**: The lower bound of the median LFC given the interval assumption (i.e., median.effect+epsilon.lower)
3. **median.effect.u**: The upper bound of the median LFC given the interval assumption (i.e., median.effect+epsilon.upper)
4. **ctt.pval**: The p-value for the composite t-test
5. **cwt.pval**: The p-value for the composite wilcoxon test
6. **gtt.pval**: The p-value for the generalized t-test
7. **ctt.pval.BH.adj**: The BH adjusted p-value for the composite t-test
8. **cwt.pval.BH.adj**: The BH adjusted p-value for the composite wilcoxon test
9. **gtt.pval.BH.adj**: The BH adjusted p-value for the generalized t-test

