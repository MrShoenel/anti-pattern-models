-   [Overview](#overview)
-   [The process model](#the-process-model)
    -   [Min/max areas between](#minmax-areas-between)
-   [The random processes](#the-random-processes)
    -   [Function to get smoothed random
        process](#function-to-get-smoothed-random-process)
    -   [Test of random process](#test-of-random-process)
    -   [Simulation of scores](#simulation-of-scores)
    -   [Examine scores](#examine-scores)
-   [Transformation of scores](#transformation-of-scores)
    -   [Empirical proof](#empirical-proof)
    -   [Visual verification](#visual-verification)

# Overview

In this notebook, we will pick an exemplary process model, and then
simulate random processes and their score with the PM, in order to learn
an empirical distribution. Then using this ePDF/eCDF, we will attempt to
properly translate and scale scores to have a uniform distribution over
the interval \[0,1\].

# The process model

We pick some linear function over the unit square as PM, and as score we
will later compute the area between the PM and each simulated random
process. The PM is shown in figure .

``` r
PM <- function(x) 0.2 + 0.2 * x

curve(expr = PM, from = 0, to = 1, ylim = c(0, 1))
```

![The exemplary process model, which is just a linear function in the
unit square.](rectify-score_files/figure-gfm/example-pm-1.png)

## Min/max areas between

From that PM we can directly observe the minimum and maximum possible
areas between. While the minimum is 0, the maximum is the area above, or
1 minus the area below. The area below is 1 × 0.2 + 1 × 0.2 ÷ 2 = 0.3:

``` r
cubature::cubintegrate(f = PM, lower = 0, upper = 1)$integral
```

    ## [1] 0.3

Therefore, the maximum distance is 0.7.

# The random processes

The distribution of the scores in an actual model once it is set up
depends on the available degrees of freedom of the model, the constant
PM itself, and the other objectives. For this experiment, we have no
such DOFs, so we will simulate processes by smoothing out some curve in
the unit squre.

## Function to get smoothed random process

Using random (but sorted) points in the unit square, we will approximate
a smooth line through these, which will become a random process.

``` r
get_smoothed_curve <- function(seed = NA, npoints = 15) {
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }

  x <- sort(c(0, 1, runif(npoints - 2)))
  y <- runif(length(x))
  temp <- loess.smooth(x = x, y = y, span = 0.35, family = "g", evaluation = 1000)
  appr <- stats::approxfun(x = ((temp$x - min(temp$x))/(max(temp$x) - min(temp$x))),
    y = temp$y/(max(temp$x) - min(temp$x)), yleft = utils::head(temp$y, 1), yright = utils::tail(temp$y,
      1))

  tempf <- Vectorize(function(x) {
    # Limit the resulting function to the bounding box of [0,1]
    min(1, max(0, appr(x)))
  })

  `attributes<-`(tempf, list(x = x, y = y))
}
```

## Test of random process

Let’s make a test, showing the PM, the random points (grey), and the
smoothed curve (red) (figure ):

``` r
test <- get_smoothed_curve(seed = 1337)

curve(PM, from = 0, to = 1, ylim = c(0, 1), lwd = 2)
lines(attr(test, "x"), attr(test, "y"), type = "l", col = "#888888")
curve2(test, from = 0, to = 1, add = TRUE, col = "red")
```

![A random process as a smooth approximation of some randomly picked
points.](rectify-score_files/figure-gfm/pm-and-observation-test-1.png)

Now the area between the two processes is:

``` r
area_diff_2_functions_score(numSamples = 1000, useYRange = c(0, 1))(PM, test)
```

    ## [1] 0.7195833

``` r
area_diff_2_functions(f1 = PM, f2 = test)
```

    ## $areas
    ## [1] 0.19100119 0.06385291 0.02556260
    ## 
    ## $value
    ## [1] 0.2804167
    ## 
    ## $intersections
    ## [1] 0.0000000 0.5938207 0.9026300 1.0000000

## Simulation of scores

Now it’s time to compute some random scores using random processes. The
more we compute, the better the approximation of the score’s real
distribution, given the degrees of freedom.

``` r
scores <- loadResultsOrCompute(file = "../results/rectify_scores.rds", computeExpr = {
  doWithParallelCluster(expr = {
    library(foreach)
    foreach::foreach(seed = seq(from = 0, to = 10000), .combine = rbind, .inorder = FALSE) %dopar%
      {
        proc = get_smoothed_curve(seed = seed)
        matrix(nrow = 1, ncol = 2, data = c(seed, area_diff_2_functions_score(numSamples = 1000,
          useYRange = c(0, 1))(f1 = PM, f2 = proc)))
      }
  })
})
```

## Examine scores

Let’s look at the distribution of scores and some other properties. The
min/max are:

``` r
scores.minmax <- c(min(scores[, 2]), max(scores[, 2]))
scores.minmax
```

    ## [1] 0.4904602 0.9379748

That means, that smallest and largest areas observed between PM and some
process were:

``` r
1 - scores.minmax
```

    ## [1] 0.50953978 0.06202516

Let’s show the best (blue) and worst (red) processes (figure ):

``` r
curve(expr = PM, from = 0, to = 1, ylim = c(0, 1), lwd = 2)
p_best <- get_smoothed_curve(seed = scores[which.max(scores[, 2]), 1])
curve2(func = p_best, from = 0, to = 1, col = "blue", add = TRUE)
p_worst <- get_smoothed_curve(seed = scores[which.min(scores[, 2]), 1])
curve2(func = p_worst, from = 0, to = 1, col = "red", add = TRUE)
```

![The best and worst processes found by our empirical
simulation.](rectify-score_files/figure-gfm/empirical-best-worst-1.png)

The empirical distributions for the scores are shown in figure :

``` r
par(mfrow = c(1, 2))

plot(stats::density(scores[, 2]))
plot(stats::ecdf(scores[, 2]))
```

![The empirical probability density and cumulative probability densities
of the simulated
scores.](rectify-score_files/figure-gfm/scores-epdf-ecdf-1.png)

The mode of the ePDF is at  ≈ 0.745 with a value of  ≈ 5.369:

``` r
temp <- stats::density(scores[, 2])
c(temp$x[which.max(temp$y)], temp$y[which.max(temp$y)])
```

    ## [1] 0.7445607 5.3686763

That means, that the relative most likely scores is  ≈ 0.745. Both, ePDF
and eCDF reveal that the distribution does not follow a uniform
distribution.

The exact quantiles are:

``` r
quantile(scores[, 2])
```

    ##        0%       25%       50%       75%      100% 
    ## 0.4904602 0.6966105 0.7480557 0.7951352 0.9379748

# Transformation of scores

We have now collected some data that indicates that the scores of our
little experiment are not uniformly distributed. Also, they do not start
at 0 and do not end at 1. The most likely score tells us that the most
likely area between both processes is  ≈ 0.255. Since the scores are not
uniformly scaled, a score just 0.1 larger than the most likely of
 ≈ 0.75 must be considered significantly better than just 0.1, since it
is much less likely:

``` r
scores.epdf <- stats::approxfun(x = temp$x, y = temp$y, yleft = 0, yright = 0)
scores.epdf(c(0.75, 0.85))
```

    ## [1] 5.340594 2.026995

Obviously, we cannot just linearly transform the scores, as it is a
non-linear relationship, see figure . However, we can simply plug in
each raw score into the eCDF to obtain a non-linear mapping to a
uniformly distributed score. In order to change the distribution of
scores to a uniform distribution *within* the observed range, one would
define the scaling function as “min(range) + ecdf(x) \* extent(range).”
However in our case, we want to obtain a score in the range \[0,1\]. So
the minimum of the target range is 0, and the extend is 1. So the
scaling will just be ecdf(x). Any eCDF is a monotonically increasing
function (often even strictly monotonically increasing), which means
that it will preserve the order between scores. So in the worst case,
two raw scores might end up closer or farther from each other, but their
order will not switch. With each additional observed value, the eCDF
will increase precision.

``` r
scores.ecdf <- stats::ecdf(scores[, 2])
curve2(func = scores.ecdf, from = 0, to = 1, lwd = 2)
tempf.lin <- function(x) x
curve2(func = tempf.lin, from = 0, to = 1, col = "red", add = TRUE)
```

![The empirical CDF and uniform CDF. The scores’ eCDF is clearly
non-linear.](rectify-score_files/figure-gfm/ecdf-vs-uniform-1.png)

The red curve is the eCDF of a uniform distribution, and the the ePDF of
our scores deviates from that significantly. First of all, there are no
scores approximately less than 0.5 or larger than roughly 0.95. The most
common score of 0.75 should be mapped to an **actual** score of 0.5. The
eCDF correctly maps it to **0.50005**. We then observe an intersection
of both CDFs at approximately 0.8. Before that, the scores’ CDF is less
than the uniform CDF, after that it is larger. That means that scores
below 0.75 are less good than they actually are (which is true), and
above that they’re better.

## Empirical proof

We can actually prove that the scores’ eCDF creates a uniformly scaled
score in range \[0,1\], by uniformly sampling from the *probability
integral transform* of the scores’ eCDF. There is no built-in way in `R`
for obtaining the inverse CDF (sometimes called “percent point function”
or PPF) of an empirical distribution, so we have to approximate our own.
For that, we’ll sample from eCDF, then invert `x` and `y`, and
approximate a function using these vectors (figure ):

``` r
temp.inv <- sapply(X = seq(min(scores[, 2]), max(scores[, 2]), length.out = 1e+05),
  FUN = scores.ecdf)
temp.inv <- (temp.inv - min(temp.inv))/(max(temp.inv) - min(temp.inv))
tempf.inv <- stats::approxfun(x = temp.inv, y = seq(min(scores[, 2]), max(scores[,
  2]), length.out = 1e+05))
curve2(tempf.inv, 0, 1)
```

![Inverse CDF (PPF) of the empirical
scores.](rectify-score_files/figure-gfm/inv-cdf-1.png)

In the following plot, we should observe an approximately uniform
distribution (figure ):

``` r
# Sampling uniformly from above PPF:
data <- tempf.inv(runif(5e+05))

par(mfrow = c(1, 2))
# First show a 're-constituted' ePDF, i.e., no non-linear mapping:
plot(stats::density(data, cut = TRUE))
# Now plug these samples into the scores' eCDF and show the distribution:
plot(stats::density(scores.ecdf(data), cut = TRUE))
```

![‘Reconstituted’ ePDF and uniform CDF of the transformed
scores.](rectify-score_files/figure-gfm/reconstructed-pdf-uniform-cdf-1.png)

Indeed, it works. Also in the left plot, we’ll see how the probability
integral transform nicely recreates a distribution very similar to the
scores’ original ePDF.

Let’s show in table where we sample from the scores’ actual range, and
how each raw score maps to a rectified score:

``` r
scores.raw <- seq(from = floor(scores.minmax[1] * 100)/100, to = ceiling(scores.minmax[2] *
  100)/100, by = 0.01)
temp <- `colnames<-`(matrix(ncol = 2, nrow = length(scores.raw), byrow = FALSE, data = c(scores.raw,
  scores.ecdf(scores.raw))), c("Raw Score", "Transformed Score"))
```

| Raw Score | Transformed Score |
|----------:|------------------:|
|      0.49 |         0.0000000 |
|      0.50 |         0.0005000 |
|      0.51 |         0.0007999 |
|      0.52 |         0.0014999 |
|      0.53 |         0.0031997 |
|      0.54 |         0.0046995 |
|      0.55 |         0.0064994 |
|      0.56 |         0.0092991 |
|      0.57 |         0.0127987 |
|      0.58 |         0.0180982 |
|      0.59 |         0.0225977 |
|      0.60 |         0.0287971 |
|      0.61 |         0.0375962 |
|      0.62 |         0.0487951 |
|      0.63 |         0.0660934 |
|      0.64 |         0.0824918 |
|      0.65 |         0.1029897 |
|      0.66 |         0.1247875 |
|      0.67 |         0.1525847 |
|      0.68 |         0.1872813 |
|      0.69 |         0.2234777 |
|      0.70 |         0.2635736 |
|      0.71 |         0.3073693 |
|      0.72 |         0.3536646 |
|      0.73 |         0.4037596 |
|      0.74 |         0.4550545 |
|      0.75 |         0.5116488 |
|      0.76 |         0.5624438 |
|      0.77 |         0.6140386 |
|      0.78 |         0.6709329 |
|      0.79 |         0.7218278 |
|      0.80 |         0.7729227 |
|      0.81 |         0.8168183 |
|      0.82 |         0.8609139 |
|      0.83 |         0.8945105 |
|      0.84 |         0.9209079 |
|      0.85 |         0.9438056 |
|      0.86 |         0.9616038 |
|      0.87 |         0.9758024 |
|      0.88 |         0.9848015 |
|      0.89 |         0.9925007 |
|      0.90 |         0.9964004 |
|      0.91 |         0.9979002 |
|      0.92 |         0.9990001 |
|      0.93 |         0.9996000 |
|      0.94 |         1.0000000 |

Raw vs. transformed scores, using the scores’ eCDF.

In the first column with the raw score, we use equally long increments
of 0.01. In the right column, the increments do not linearly correspond
to that.

## Visual verification

Let’s generate some more random processes and show their raw and
calibrated area-in-between scores (figure ).

``` r
par(mfrow = c(3, 3))

for (i in 1:9) {
  P <- get_smoothed_curve(seed = 6e+05 + i)
  score_raw <- area_diff_2_functions_score(useYRange = c(0, 1))(f1 = PM, f2 = P)
  score_cal <- scores.ecdf(score_raw)
  curve2(func = PM, from = 0, to = 1, lwd = 2, ylim = c(0, 1), xlab = "", ylab = "",
    main = paste0(i), sub = paste0("raw=", format(score_raw, nsmall = 3, digits = 3),
      "; ", "cal=", format(score_cal, nsmall = 3, digits = 3), "; ", "area=",
      format(1 - score_raw, nsmall = 3, digits = 3)))
  curve2(func = P, from = 0, to = 1, col = "red", add = TRUE)
}
```

![Some random processes, their raw and calibrated scores. Also shown is
the area in
between.](rectify-score_files/figure-gfm/visual-verify-1.png)

The random processes of figure reflect almost the entire spectrum of all
previously observed raw scores, we almost see the lowest lows and
highest highs. Recall that the score computed is 1 minus the area in
between processes (which can maximally be 0.7). Looking at examples 1,
5, 6, and 7 we observe a rather large area, and while the raw scores do
not look so bad, the calibrated scores reveal the bad fit of PM/P.
