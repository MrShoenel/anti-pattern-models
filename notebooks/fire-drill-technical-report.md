





-   [Introduction](#introduction)
-   [Data](#data)
    -   [The Ground Truth](#the-ground-truth)
    -   [The Student Projects](#the-student-projects)
    -   [Modeling of metrics and events using
        KDE](#modeling-of-metrics-and-events-using-kde)
-   [Patterns for scoring the
    projects](#patterns-for-scoring-the-projects)
    -   [Pattern I: Initial best guess](#pattern-i-initial-best-guess)
        -   [Initialize the pattern](#initialize-the-pattern)
    -   [Pattern II: Adaptation of best
        guess](#pattern-ii-adaptation-of-best-guess)
        -   [Preparation](#preparation)
        -   [Defining the losses](#defining-the-losses)
        -   [Fitting the pattern](#fitting-the-pattern)
        -   [Inversing the parameters](#inversing-the-parameters)
    -   [Pattern III: Averaging the ground
        truth](#pattern-iii-averaging-the-ground-truth)
    -   [Pattern III (b): Evidence-based](#pattern-iii-b-evidence-based)
        -   [Preparation](#preparation-1)
        -   [Finding the best fit](#finding-the-best-fit)
        -   [Create pattern from best
            fit](#create-pattern-from-best-fit)
-   [Scoring of projects](#scoring-of-projects)
    -   [Scoring mechanisms](#scoring-mechanisms)
    -   [Pattern I](#pattern-i)
    -   [Pattern II](#pattern-ii)
    -   [Pattern II (without alignment)](#pattern-ii-without-alignment)
    -   [Pattern III (average)](#pattern-iii-average)
    -   [Pattern III (average, no
        alignment)](#pattern-iii-average-no-alignment)
        -   [Linear combination of
            scores](#linear-combination-of-scores)
    -   [Pattern III (b)](#pattern-iii-b)
        -   [Linear combination of
            scores](#linear-combination-of-scores-1)
-   [References](#references)

Introduction
============

This is the complementary technical report for the paper/article
tentatively entitled “Towards Data-Based Detection of Project Management
Anti-Patterns”. Here, we import the ground truth as well as all
projects’ data, and instantiate our model based on *self-regularizing
Boundary Time Warping and Boundary Amplitude Warping*. Given a few
patterns that represent the **Fire Drill** anti-pattern (AP), the goal
is evaluate these patterns and their aptitude for detecting the AP in
concordance with the ground truth.

All complementary data and results can be found at Zenodo. This notebook
was written in a way that it can be run without any additional efforts
to reproduce the outputs (using the pre-computed results). This notebook
has a canonical
URL<sup>[\[Link\]](https://github.com/sse-lnu/anti-pattern-models/blob/master/notebooks/fire-drill-technical-report.Rmd)</sup>
and can be read online as a rendered
markdown<sup>[\[Link\]](https://github.com/sse-lnu/anti-pattern-models/blob/master/notebooks/fire-drill-technical-report.md)</sup>
version. All code can be found in this repository, too.

Data
====

We have 9 projects conducted by students, and two raters have
**independently**, i.e., without prior communication, assessed to what
degree the AP is present in each project. This was done using a scale
from zero to ten, where zero means that the AP was not present, and ten
would indicate a strong manifestation.

The Ground Truth
----------------

``` r
ground_truth <- read.csv(file = "../data/ground-truth.csv", sep = ";")
```

| project    |  rater.a|  rater.b|  consensus|  rater.mean|
|:-----------|--------:|--------:|----------:|-----------:|
| project\_1 |        2|        0|          1|         1.0|
| project\_2 |        0|        0|          0|         0.0|
| project\_3 |        8|        5|          6|         6.5|
| project\_4 |        8|        6|          8|         7.0|
| project\_5 |        1|        1|          1|         1.0|
| project\_6 |        4|        1|          2|         2.5|
| project\_7 |        2|        3|          3|         2.5|
| project\_8 |        0|        0|          0|         0.0|
| project\_9 |        1|        4|          5|         2.5|

Using the *quadratic weighted Kappa* (Cohen 1968), we can report an
unadjusted agreement of **0.715** for both raters. A Kappa value in the
range \[0.6, 0.8\] is considered *substantial*, and values beyond that
as *almost perfect* (Landis and Koch 1977). As for the
Pearson-correlation, we report a slightly higher value of **0.771**. The
entire ground truth is shown in table . The final consensus was reached
after both raters exchanged their opinions, and it is the consensus that
we will use as the actual ground truth from here on and out.

The Student Projects
--------------------

The ground truth was extracted from nine student-conducted projects.
Seven of these were implemented simultaneously between March and June
2020, and two the year before in a similar timeframe.

``` r
student_projects <- read.csv(file = "../data/student-projects.csv", sep = ";")
```

We have a total of:

-   Nine projects,
-   37 authors that authored 1219 commits total which are of type
-   Adaptive / Corrective / Perfective (`a/c/p`) commits: 392 / 416 /
    411

We have a complete breakdown of all activities across all projects in
figure .

``` r
student_projects_info <- NULL

for (pId in unique(student_projects$project)) {
  temp <- student_projects[student_projects$project == pId, ]
  student_projects_info <- rbind(student_projects_info, data.frame(
    project = pId,
    authors = length(unique(temp$AuthorNominalLabel)),
    commits = nrow(temp),
    a = nrow(temp[temp$label == "a", ]),
    c = nrow(temp[temp$label == "c", ]),
    p = nrow(temp[temp$label == "p", ]),
    avgDens = round(mean(temp$Density), 3)
  ))
}
```

| project    |  authors|  commits|    a|    c|    p|  avgDens|
|:-----------|--------:|--------:|----:|----:|----:|--------:|
| project\_1 |        4|      116|   36|   32|   48|    0.879|
| project\_2 |        5|      226|   42|  108|   76|    0.891|
| project\_3 |        4|      111|   26|   35|   50|    0.785|
| project\_4 |        4|      126|   29|   59|   38|    0.870|
| project\_5 |        2|      110|   33|   26|   51|    0.814|
| project\_6 |        4|      217|   79|   63|   75|    0.784|
| project\_7 |        5|      183|   83|   64|   36|    0.813|
| project\_8 |        4|       30|   10|    6|   14|    0.687|
| project\_9 |        5|      100|   54|   23|   23|    0.743|

<img src="fire-drill-technical-report_files/figure-markdown_github/project-activity-1.png" alt="Commit activities across projects"  />
<p class="caption">
Commit activities across projects
</p>

We have slightly different begin- and end-times in each project.
However, the data for all projects was previously cropped, so that each
project’s extent marks the absolute begin and end of it – it starts with
the first commit and ends with the last. As for our methods here, we
only need to make sure that we scale the timestamps into a relative
\[0, 1\]-range, where 1 marks the project’s end.

For each project, we model **four** variables: The activities
**adaptive** (**`A`**), **corrective+perfective** (**`CP`**), the
frequency of all activities, regardless of their type (**`FREQ`**), and
the **source code density** (**`SCD`**). While for the first three
variables we estimate a Kernel density, the last variable is a metric
collected with each commit. The data for it is mined using `Git-Density`
(Hönel 2020), and we use a highly efficient commit classification
model[1] ( ≈ 83.6% accuracy,  ≈ 0.745 Kappa) (Hönel et al. 2020) to
attach maintenance activity labels to each commit, based on size- and
keyword-data only.

Technically, we will compose each variable into an instance of our
`Signal`-class. Before we start, we will do some normalizations and
conversions, like converting the timestamps. This has to be done on a
per-project basis.

``` r
student_projects$label <- as.factor(student_projects$label)
student_projects$project <- as.factor(student_projects$project)
student_projects$AuthorTimeNormalized <- NA_real_
```

``` r
for (pId in levels(student_projects$project)) {
  student_projects[student_projects$project == pId, ]$AuthorTimeNormalized <-
    (student_projects[student_projects$project == pId, ]$AuthorTimeUnixEpochMilliSecs -
      min(student_projects[student_projects$project == pId, ]$AuthorTimeUnixEpochMilliSecs))
  student_projects[student_projects$project == pId, ]$AuthorTimeNormalized <-
    (student_projects[student_projects$project == pId, ]$AuthorTimeNormalized /
      max(student_projects[student_projects$project == pId, ]$AuthorTimeNormalized))
}
```

And now for the actual signals: Since the timestamps have been
normalized for each project, we model each variable to actually start at
0 and end at 1 (the support). We will begin with activity-related
variables before we model the source code density, as the process is
different. When using Kernel density estimation (KDE), we obtain an
empirical probability density function (PDF) that integrates to 1. This
is fine when looking at all activities combined (**`FREQ`**). However,
when we are interested in a specific fraction of the activities, say
**`A`**, then we should scale its activities according to its overall
ratio. Adding all scaled activities together should again integrate to
1. When this is done, we scale one last time such that no empirical PDF
has a co-domain larger than 1.

``` r
project_signals <- list()

# passed to stats::density
use_kernel <- "gauss" # "rect"

for (pId in levels(student_projects$project)) {
  temp <- student_projects[student_projects$project == pId, ]
  
  # We'll need these for the densities:
  acp_ratios <- table(temp$label) / sum(table(temp$label))
  
  dens_a <- densitySafe(
    from = 0, to = 1, safeVal = NA_real_,
    data = temp[temp$label == "a", ]$AuthorTimeNormalized,
    ratio = acp_ratios[["a"]], kernel = use_kernel)
  
  dens_cp <- densitySafe(
    from = 0, to = 1, safeVal = NA_real_,
    data = temp[temp$label == "c" | temp$label == "p", ]$AuthorTimeNormalized,
    ratio = acp_ratios[["c"]] + acp_ratios[["p"]], kernel = use_kernel)
  
  dens_freq <- densitySafe(
    from = 0, to = 1, safeVal = NA_real_,
    data = temp$AuthorTimeNormalized, ratio = 1, kernel = use_kernel)
  
  # All densities need to be scaled together once more, by dividing
  # for the maximum value of the FREQ-variable.
  ymax <- max(c(attr(dens_a, "ymax"), attr(dens_cp, "ymax"), attr(dens_freq, "ymax")))
  dens_a <- stats::approxfun(
    x = attr(dens_a, "x"), y = sapply(X = attr(dens_a, "x"), FUN = dens_a) / ymax)
  dens_cp <- stats::approxfun(
    x = attr(dens_cp, "x"), y = sapply(X = attr(dens_cp, "x"), FUN = dens_cp) / ymax)
  dens_freq <- stats::approxfun(
    x = attr(dens_freq, "x"), y = sapply(X = attr(dens_freq, "x"), FUN = dens_freq) / ymax)
  
  project_signals[[pId]] <- list(
    A = Signal$new(name = paste(pId, "A", sep = "_"),
                   func = dens_a, support = c(0, 1), isWp = FALSE),
    CP = Signal$new(name = paste(pId, "CP", sep = "_"),
                    func = dens_cp, support = c(0, 1), isWp = FALSE),
    FREQ = Signal$new(name = paste(pId, "FREQ", sep = "_"),
                      func = dens_freq, support = c(0, 1), isWp = FALSE)
  )
}
```

Now, for each project, we estimate the variable for the source code
density as follows:

``` r
for (pId in levels(student_projects$project)) {
  temp <- data.frame(
    x = student_projects[student_projects$project == pId, ]$AuthorTimeNormalized,
    y = student_projects[student_projects$project == pId, ]$Density)
  temp <- temp[with(temp, order(x)), ]
  
  # Using a polynomial with maximum possible degree, we smooth the
  # SCD-data, as it can be quite "peaky"
  temp_poly <- poly_autofit_max(x = temp$x, y = temp$y, startDeg = 13)
  
  dens_scd <- Vectorize((function() {
    rx <- range(temp$x)
    ry <- range(temp$y)
    poly_y <- stats::predict(temp_poly, x = temp$x)
    tempf <- stats::approxfun(x = temp$x, y = poly_y, ties = "ordered")
    function(x) {
      if (x < rx[1] || x > rx[2]) {
        return(NA_real_)
      }
      max(ry[1], min(ry[2], tempf(x)))
    }
  })())
  
  project_signals[[pId]][["SCD"]] <- Signal$new(
    name = paste(pId, "SCD", sep = "_"), func = dens_scd,
    support = c(0, 1), isWp = FALSE)
}
```

Let’s plot all the projects:

<img src="fire-drill-technical-report_files/figure-markdown_github/project-vars-1.png" alt="All variables over each project's time span"  />
<p class="caption">
All variables over each project’s time span
</p>

Modeling of metrics and events using KDE
----------------------------------------

We need to make an important distinction between events and metrics. An
event does not carry other information, other than that it occurred. One
could thus say that such an event is *nulli*-variate. If an event were
to carry extra information, such as a measurement that was taken, it
would be *uni*-variate. That is the case for many metrics in software:
the time of their measurement coincides with an event, such as a commit
that was made. On the time-axis we thus know **when** it occurred and
**what** was its value. Such a metric could be easily understood as a
*bivariate x/y* variable and be plotted in a two-dimensional Cartesian
coordinate system.

An event however does not have any associated y-value we could plot.
Given a time-axis, we could make a mark whenever it occurred. Some of
the markers would probably be closer to each other or be more or less
accumulated. The y-value could express these accumulations relative to
each other. These are called *densities*. This is exactly what KDE does:
it expresses the relative accumulations of data on the x-axis as density
on the y-axis. For KDE, the actual values on the x-axis have another
meaning, and that is to compare the relative likelihoods of the values
on it, since the axis is ordered. For our case however, the axis is
linear time and carries no such meaning. The project data we analyze is
a kind of sampling over the project’s events. We subdivide the gathered
project data hence into these two types of data series:

-   **Events**: They do not carry any extra information or measurements.
    As for the projects we analyze, events usually are occurrences of
    specific types of commits, for example. The time of occurrence is
    the x-value on the time-axis, and the y-value is obtained through
    KDE. We model all maintenance activities as such variables.
-   **Metrics**: Extracted from the project at specific times, for
    example at every commit. We can extract any number or type of
    metric, but each becomes its own variable, where the x-value is on
    the time-axis, and the y-value is the metric’s value. We model the
    source code density as such a variable.

Patterns for scoring the projects
=================================

Our overall goal is to propose a single model that is able to detect the
presence of the Fire Drill AP, and how strong its manifestation is. In
order to do that, we require a pattern that defines how a Fire Drill
looks in practice. Any real-world project can never follow such a
pattern perfectly, because of, e.g., time dilation and compression. Even
after correcting these, some distance between the project and the
pattern will remain. The projects from figure indicate that certain
phases occur, but that their occurrence happens at different points in
time, and lasts for various durations.

Given some pattern, we first attempt to remove any distortions in the
data, by using our new model *self-regularizing Boundary Time Warping*
(sr-BTW). This model takes a pattern that is subdivided into one or more
intervals, and aligns the project data such that the loss in each
interval is minimized. After alignment, we calculate a score that
quantifies the remaining differences. Ideally, we hope to find a
(strong) positive correlation of these scores with the ground truth.

Pattern I: Initial best guess
-----------------------------

``` r
fd_data_concat <- readRDS("../data/fd_data_concat.rds")
```

This pattern was created based on all available literature, **without**
inspecting any of the projects. It is subdivided into four intervals:

1.  Begin – Short project warm-up phase
2.  Long Stretch – The longest phase in the project, about which we do
    not know much about, except for that there should be a rather
    constant amount of activities over time.
3.  Fire Drill – Characteristic is a sudden and steep increase of
    adaptive activities. This phase is over once these activities
    reached their apex.
4.  Aftermath – Everything after the apex. We should see even steeper
    declines.

Brown et al. (1998) describe a typical scenario where about six months
are spent on non-developmental activities, and the actual software is
then developed in less than four weeks. If we were to include some of
the aftermath, the above first guess would describe a project of about
eight weeks.

We define the boundaries as follows (there are three boundaries to split
the pattern into four intervals):

``` r
fd_data_boundaries <- c("b1" = 0.085, "b2" = 0.625, "b3" = 0.875)
```

The pattern and its boundaries look like this:

``` r
plot_project_data(data = fd_data_concat, boundaries = fd_data_boundaries)
```

<img src="fire-drill-technical-report_files/figure-markdown_github/pattern-1-1.png" alt="The pattern that was our initial best guess"  />
<p class="caption">
The pattern that was our initial best guess
</p>

### Initialize the pattern

The pattern as shown in is just a collection of x/y coordinate-data, and
for us being able to use it, we need to instantiate it. We do this by
storing each variable in an instance of `Signal`.

``` r
p1_signals <- list(
  A = Signal$new(name = "p1_A", support = c(0, 1), isWp = TRUE, func = 
    stats::approxfun(
      x = fd_data_concat[fd_data_concat$t == "A", ]$x,
      y = fd_data_concat[fd_data_concat$t == "A", ]$y)),
  CP = Signal$new(name = "p1_CP", support = c(0, 1), isWp = TRUE, func = 
    stats::approxfun(
      x = fd_data_concat[fd_data_concat$t == "CP", ]$x,
      y = fd_data_concat[fd_data_concat$t == "CP", ]$y)),
  FREQ = Signal$new(name = "p1_FREQ", support = c(0, 1), isWp = TRUE, func = 
    stats::approxfun(
      x = fd_data_concat[fd_data_concat$t == "FREQ", ]$x,
      y = fd_data_concat[fd_data_concat$t == "FREQ", ]$y)),
  SCD = Signal$new(name = "p1_SCD", support = c(0, 1), isWp = TRUE, func = 
    stats::approxfun(
      x = fd_data_concat[fd_data_concat$t == "SCD", ]$x,
      y = fd_data_concat[fd_data_concat$t == "SCD", ]$y))
)
```

<img src="fire-drill-technical-report_files/figure-markdown_github/unnamed-chunk-17-1.png" alt="The separate signals of pattern I."  />
<p class="caption">
The separate signals of pattern I.
</p>

Pattern II: Adaptation of best guess
------------------------------------

The second pattern is a compromise between the first and the third:
While we want to keep as much of the initial best guess, we also want to
adjust the pattern based on the projects and the ground truth. Adjusting
means, that we will keep what is in each interval, but we allow each
interval to stretch and compress, and we allow each interval to impose a
vertical translation both at then begin and end (a somewhat trapezoidal
translation). In any case, each such alteration is a linear affine
transformation. Additionally to sr-BTW, we will also apply **sr-BAW**
(self-regularizing Boundary Amplitude Warping) to accomplish this. This
model is called **`srBTAW`** and the process is the following:

-   The pattern is decomposed into its four variables first, as we can
    adapt these (almost) independently from each other.
-   Then, for each type of variable, an instance of `srBTAW` is created.
    As **Warping Candidates** (WC) we add all of the projects’
    corresponding variables. The **Warping Pattern** (WP) is the single
    variable from the pattern in this case – again, we warp the project
    data, however, eventually the learned warping gets inversed and
    applied to the WC.
-   All four `srBTAW` instances are then fitted simultaneously: While we
    allow the y-translations to adapt independently for each type of
    variable, all instances share the same intervals, as eventually we
    have to assemble the variables back into a common pattern.

### Preparation

We already have the `srBTAW` **Multilevel model**, which can keep track
of arbitrary many variables and losses. The intention behind this
however was, to track variables of the **same type**, i.e., signals that
are logically of the same type. In our case this means that any single
instance should only track variables that are either `A`, `CP`, `FREQ`
or `SCD`. For this pattern, the WP is a single signal per variable, and
the WC is the corresponding signal from each of the nine projects. This
is furthermore important to give different weights to different
variables. In our case, we want to give a lower weight to the
`SCD`-variable.

As for the loss, we will first test a combined loss that measures
**`3`** properties: The area between curves (or alternatively the
residual sum of squares), the correlation between the curves, and the
arc-length ratio between the curves. We will consider any of these to be
equally important, i.e., no additional weights. Each loss shall cover
all intervals with weight  = 1, except for the Long Stretch interval,
where we will use a reduced weight.

There are 4 types of variables, 7 projects (two projects have consensus
 = 0, i.e., no weight) and 2 × 3 single losses, resulting in 168 losses
to compute. The final weight for each loss is computed as:
*ω*<sub>*i*</sub> = *ω*<sup>(project)</sup> × *ω*<sup>(vartype)</sup> × *ω*<sup>(interval)</sup>.
For the phase Long Stretch, the weight for any loss will $\\frac{1}{2}$,
and for the source code density we will chose $\\frac{1}{2}$, too. The
weight of each project is based on the consensus of the ground truth.
The ordinal scale for that is \[0, 10\], so that we will divide the
score by 10 and use that as weight. Examples:

-   **A** in Fire Drill in project *p*3: *ω* = 0.6 × 1 × 1 = 0.6
    (consensus is 6 in project *p*3)
-   **FREQ** in Long Stretch in project *p*7: *ω* = 0.3 × 0.5 × 1 = 0.15
    and
-   **SCD** in Long Stretch in project *p*4:
    *ω* = 0.8 × 0.5 × 0.5 = 0.2.

In table we show all projects with a consensus-score  \> 0, projects 2
and 8 are not included any longer.

``` r
ground_truth$consensus_score <- ground_truth$consensus / 10
weight_vartype <- c("A" = 1, "CP" = 1, "FREQ" = 1, "SCD" = 0.5)
weight_interval <- c("Begin" = 1, "Long Stretch" = 0.5, "Fire Drill" = 1, "Aftermath" = 1)
```

``` r
temp <- expand.grid(weight_interval, weight_vartype, ground_truth$consensus_score)
temp$p <- temp$Var1 * temp$Var2 * temp$Var3
weight_total <- sum(temp$p)
```

The sum of all weights combined is 31.85.

|     | project    |  consensus|  consensus\_score|
|:----|:-----------|----------:|-----------------:|
| 1   | project\_1 |          1|               0.1|
| 3   | project\_3 |          6|               0.6|
| 4   | project\_4 |          8|               0.8|
| 5   | project\_5 |          1|               0.1|
| 6   | project\_6 |          2|               0.2|
| 7   | project\_7 |          3|               0.3|
| 9   | project\_9 |          5|               0.5|

### Defining the losses

For the optimization we will use mainly **`5`** classes:

-   `srBTAW_MultiVartype`: One instance globally, that manages all
    parameters across all instances of `srBTAW`.
-   `srBTAW`: One instance per variable-type, so here we’ll end up with
    four instances.
-   `srBTAW_LossLinearScalarizer`: A linear scalarizer that will take on
    all of the defined singular losses and compute and add them together
    according to their weight.
-   `srBTAW_Loss2Curves`: Used for each of the 168 singular losses, and
    configured using a specific loss function, weight, and set of
    intervals where it ought to be used.
-   `TimeWarpRegularization`: One global instance for all `srBTAW`
    instances, to regularize extreme intervals. We chose a mild weight
    for this of just 1, which is small compared to the sum of all other
    weights (31.85).

``` r
p2_smv <- srBTAW_MultiVartype$new()

p2_vars <- c("A", "CP", "FREQ", "SCD")
p2_inst <- list()
for (name in p2_vars) {
  p2_inst[[name]] <- srBTAW$new(
    theta_b = c(0, fd_data_boundaries, 1),
    gamma_bed = c(0, 1, sqrt(.Machine$double.eps)),
    lambda = rep(sqrt(.Machine$double.eps), length(p2_vars)),
    begin = 0, end = 1, openBegin = FALSE, openEnd = FALSE,
    useAmplitudeWarping = TRUE,
    # We allow these to be larger; however, the final result should be within [0,1]
    lambda_ymin = rep(-10, length(p2_vars)),
    lambda_ymax = rep( 10, length(p2_vars)),
    isObjectiveLogarithmic = TRUE,
    paramNames = c("v",
                   paste0("vtl_", seq_len(length.out = length(p2_vars))),
                   paste0("vty_", seq_len(length.out = length(p2_vars)))))
  
  # We can already add the WP:
  p2_inst[[name]]$setSignal(signal = p1_signals[[name]])
  p2_smv$setSrbtaw(varName = name, srbtaw = p2_inst[[name]])
  
  # .. and also all the projects' signals:
  for (project in ground_truth[ground_truth$consensus > 0, ]$project) {
    p2_inst[[name]]$setSignal(signal = project_signals[[project]][[name]])
  }
}

# We call this there so there are parameters present.
set.seed(1337)
p2_smv$setParams(params =
  `names<-`(x = runif(n = p2_smv$getNumParams()), value = p2_smv$getParamNames()))
```

We can already initialize the linear scalarizer. This includes also to
set up some progress-callback. Even with massive parallelization, this
process will take its time so it will be good to know where we are
approximately.

``` r
p2_lls <- srBTAW_LossLinearScalarizer$new(
  returnRaw = FALSE,
  computeParallel = TRUE, progressCallback = function(what, step, total) {
    # if (step == total) {
    #   print(paste(what, step, total))
    # }
  })

for (name in names(p2_inst)) {
  p2_inst[[name]]$setObjective(obj = p2_lls)
}
```

The basic infrastructure stands, so now it’s time to instantiate all of
the singular losses. First we define a helper-function to do the
bulk-work, then we iterate all projects, variables and intervals.

``` r
#' This function creates a singular loss that is a linear combination
#' of an area-, correlation- and arclength-loss (all with same weight).
p2_attach_combined_loss <- function(project, vartype, intervals) {
  weight_p <- ground_truth[ground_truth$project == project, ]$consensus_score
  weight_v <- weight_vartype[[vartype]]
  temp <- weight_interval[intervals]
  stopifnot(length(unique(temp)) == 1)
  weight_i <- unique(temp)
  weight <- weight_p * weight_v * weight_i
  
  lossRss <- srBTAW_Loss_Rss$new(
    wpName = paste0("p1_", vartype), wcName = paste(project, vartype, sep = "_"),
    weight = weight, intervals = intervals, continuous = FALSE,
    numSamples = rep(500, length(intervals)), returnRaw = TRUE)
  
  p2_inst[[vartype]]$addLoss(loss = lossRss)
  p2_lls$setObjective(
    name = paste(project, vartype, paste(intervals, collapse = "_"),
                 "rss", sep = "_"), obj = lossRss)
}
```

Let’s call our helper iteratively:

``` r
interval_types <- list(A = c(1,3,4), B = 2)

for (vartype in p2_vars) {
  for (project in ground_truth[ground_truth$consensus > 0, ]$project) {
    for (intervals in interval_types) {
      p2_attach_combined_loss(
        project = project, vartype = vartype, intervals = intervals)
    }
  }
  
  # Add one per variable-type:
  lossYtrans <- YTransRegularization$new(
    wpName = paste0("p1_", vartype), wcName = paste(project, vartype, sep = "_"),
    intervals = seq_len(length.out = 4), returnRaw = TRUE,
    weight = 1, use = "tikhonov")

  p2_inst[[vartype]]$addLoss(loss = lossYtrans)
  p2_lls$setObjective(
    name = paste(vartype, "p2_reg_output", sep = "_"),
    obj = lossYtrans)
}
```

Finally, we add the regularizer for extreme intervals:

``` r
p2_lls$setObjective(name = "p2_reg_exint2", obj = TimeWarpRegularization$new(
  weight = 0.25 * p2_lls$getNumObjectives(), use = "exint2", returnRaw = TRUE,
  wpName = p1_signals$A$getName(), wcName = project_signals$project_1$A$getName(),
  intervals = seq_len(length.out = length(p2_vars))
)$setSrBtaw(srbtaw = p2_inst$A))
```

### Fitting the pattern

``` r
p2_params <- loadResultsOrCompute(file = "../results/p2_params.rds", computeExpr = {
  cl <- parallel::makePSOCKcluster(min(64, parallel::detectCores()))
  tempf <- tempfile()
  saveRDS(object = list(a = p2_smv, b = p2_lls), file = tempf)
  parallel::clusterExport(cl, varlist = list("tempf"))
  
  res <- doWithParallelClusterExplicit(cl = cl, expr = {
    optimParallel::optimParallel(
      par = p2_smv$getParams(),
      method = "L-BFGS-B",
      lower = c(
        rep(-.Machine$double.xmax, length(p2_vars)), # v_[vartype]
        rep(sqrt(.Machine$double.eps), length(p2_vars)), # vtl
        rep(-.Machine$double.xmax, length(p2_vars) * length(weight_vartype))), # vty for each
      upper = c(
        rep(.Machine$double.xmax, length(p2_vars)),
        rep(1, length(p2_vars)),
        rep(.Machine$double.xmax, length(p2_vars) * length(weight_vartype))),
      fn = function(x) {
        temp <- readRDS(file = tempf)
        temp$a$setParams(params = x)
        temp$b$compute0()
      },
      parallel = list(cl = cl, forward = FALSE, loginfo = TRUE)
    )
  })
})

p2_fr <- FitResult$new("a")
p2_fr$fromOptimParallel(p2_params)
format(p2_fr$getBest(paramName = "loss")[
  1, !(p2_fr$getParamNames() %in% c("duration", "begin", "end"))],
  scientific = FALSE, digits = 4, nsmall = 4)
```

    ##         v_A        v_CP      v_FREQ       v_SCD       vtl_1       vtl_2 
    ## "-0.014674" "-0.001527" "-0.008464" "-0.325218" " 0.486963" " 0.040402" 
    ##       vtl_3       vtl_4     vty_1_A    vty_1_CP  vty_1_FREQ   vty_1_SCD 
    ## " 0.493052" " 0.988596" " 0.055662" " 0.285855" " 0.413141" " 0.398693" 
    ##     vty_2_A    vty_2_CP  vty_2_FREQ   vty_2_SCD     vty_3_A    vty_3_CP 
    ## "-0.004618" " 0.023071" "-0.064090" "-0.149528" " 0.661491" "-0.333191" 
    ##  vty_3_FREQ   vty_3_SCD     vty_4_A    vty_4_CP  vty_4_FREQ   vty_4_SCD 
    ## " 0.462005" " 0.263016" "-1.048072" "-0.287820" "-1.437350" "-0.342797" 
    ##        loss 
    ## " 6.167237"

<img src="fire-drill-technical-report_files/figure-markdown_github/p2-params-fig-1.png" alt="Neg. Log-loss of fitting pattern type II."  />
<p class="caption">
Neg. Log-loss of fitting pattern type II.
</p>

### Inversing the parameters

For this pattern, we have warped all the projects to the pattern, while
the ultimate goal is to warp the pattern to all the projects (or,
better, to warp each type of variable of the WP to the group of
variables of the same type of all projects, according to their weight,
which is determined by the consensus of the ground truth). So, if we
know how to go from A to B, we can inverse the learned parameters and go
from B to A, which means in our case that we have to apply the inverse
parameters to the WP in order to obtain WP-prime.

As for y-translations (that is, *v*, as well as all
**ϑ**<sup>(*y*)</sup>), the inversion is simple: we multiply these
parameters with  − 1. The explanation for that is straightforward: If,
for example, we had to go down by  − 0.5, to bring the data closer to
the pattern, then that means that we have to lift the pattern by  + 0.5
to achieve the inverse effect.

Inversing the the boundaries is simple, too, and is explained by how we
take some portion of the WC (the source) and warp it to the
corresponding interval of the WP (the target).

That’s how we do it:

-   Given are the WP’s **original** boundaries, **θ**<sup>(*b*)</sup>,
    and the learned **ϑ**<sup>(*l*)</sup>. The goal is, for each *q*-th
    interval, to take what is in the WP’s interval and warp it according
    to the learned length.
-   Given the boundaries-to-lengths operator, T<sup>(*l*)</sup>, and the
    lengths-to-boundaries operator, T<sup>(*b*)</sup>, we can convert
    between **θ** and **ϑ**.
-   Start with a new instance of `SRBTW` (or `SRBTWBAW` for also warping
    y-translations) and set as
    **θ**<sup>(*b*)</sup> = T<sup>(*b*)</sup>(**ϑ**<sup>(*l*)</sup>).
    The learned lengths will become the **target** intervals.
-   Add the variable that ought to be transformed as **WC**, and set
    **ϑ**<sup>(*l*)</sup> = T<sup>(*l*)</sup>(**θ**<sup>(*b*)</sup>).
-   That will result in that we are taking what was *originally* in each
    interval, and warp it to a new length.
-   The warped signal is then the `M`-function of the
    `SRBTW`/`SRBTWBAW`-instance.

Short example: Let’s take the `SCD`-variable from the first pattern and
warp it!

``` r
# Transforming some learned lengths to new boundaries:
p2_ex_thetaB <- c(0, .3, .5, .7, 1)
# Transforming the original boundaries to lengths:
p2_ex_varthetaL <- unname(c(
  fd_data_boundaries[1], fd_data_boundaries[2] - fd_data_boundaries[1],
  fd_data_boundaries[3] - fd_data_boundaries[2], 1 - fd_data_boundaries[3]))

p2_ex_srbtw <- SRBTW$new(
  theta_b = p2_ex_thetaB, gamma_bed = c(0, 1, 0),
  wp = p1_signals$SCD$get0Function(), wc = p1_signals$SCD$get0Function(),
  lambda = rep(0, 4), begin = 0, end = 1)

p2_ex_srbtw$setParams(vartheta_l = p2_ex_varthetaL)
```

In figure we can quite clearly see how the pattern warped from the blue
intervals into the orange intervals

<img src="fire-drill-technical-report_files/figure-markdown_github/inverse-example-fig-1.png" alt="Warping the variable from within the blue to the orange intervals."  />
<p class="caption">
Warping the variable from within the blue to the orange intervals.
</p>

We have learned the following parameters from our optimization for
pattern II:

``` r
p2_best <- p2_fr$getBest(paramName = "loss")[
  1, !(p2_fr$getParamNames() %in% c("begin", "end", "duration"))]
p2_best
```

    ##          v_A         v_CP       v_FREQ        v_SCD        vtl_1        vtl_2 
    ## -0.014673620 -0.001527483 -0.008464488 -0.325217649  0.486962785  0.040401601 
    ##        vtl_3        vtl_4      vty_1_A     vty_1_CP   vty_1_FREQ    vty_1_SCD 
    ##  0.493052472  0.988595727  0.055661794  0.285855275  0.413141394  0.398693496 
    ##      vty_2_A     vty_2_CP   vty_2_FREQ    vty_2_SCD      vty_3_A     vty_3_CP 
    ## -0.004618406  0.023071288 -0.064090457 -0.149527947  0.661491448 -0.333190639 
    ##   vty_3_FREQ    vty_3_SCD      vty_4_A     vty_4_CP   vty_4_FREQ    vty_4_SCD 
    ##  0.462004516  0.263016381 -1.048071875 -0.287819945 -1.437350319 -0.342797057 
    ##         loss 
    ##  6.167237270

All of the initial translations (*v*) are zero. The learned lengths
converted to boundaries are:

``` r
# Here, we transform the learned lengths to boundaries.
p2_best_varthetaL <- p2_best[names(p2_best) %in% paste0("vtl_", 1:4)] /
  sum(p2_best[names(p2_best) %in% paste0("vtl_", 1:4)])
p2_best_varthetaL
```

    ##      vtl_1      vtl_2      vtl_3      vtl_4 
    ## 0.24238912 0.02011018 0.24542030 0.49208041

``` r
p2_best_thetaB <- unname(c(
  0, p2_best_varthetaL[1], sum(p2_best_varthetaL[1:2]),
  sum(p2_best_varthetaL[1:3]), 1))
p2_best_thetaB
```

    ## [1] 0.0000000 0.2423891 0.2624993 0.5079196 1.0000000

The first two intervals are rather short, while the last two are
comparatively long. Let’s transform all of the pattern’s variables
according to the parameters:

``` r
p2_signals <- list()

for (vartype in names(weight_vartype)) {
  temp <- SRBTWBAW$new(
    theta_b = unname(p2_best_thetaB), gamma_bed = c(0, 1, 0),
    wp = p1_signals[[vartype]]$get0Function(), wc = p1_signals[[vartype]]$get0Function(),
    lambda = rep(0, 4), begin = 0, end = 1,
    lambda_ymin = rep(0, 4), lambda_ymax = rep(1, 4)) # not important here
  # That's still the same ('p2_ex_varthetaL' is the original
  # boundaries of Pattern I transformed to lengths):
  temp$setParams(vartheta_l = p2_ex_varthetaL,
                 v = -1 * p2_best[paste0("v_", vartype)],
                 vartheta_y = -1 * p2_best[paste0("vty_", 1:4, "_", vartype)])
  
  p2_signals[[vartype]] <- Signal$new(
    name = paste0("p2_", vartype), support = c(0, 1), isWp = TRUE, func = Vectorize(temp$M))
}
```

The 2nd pattern, as derived from the ground truth, is shown in figure .

<img src="fire-drill-technical-report_files/figure-markdown_github/p2-signals-1.png" alt="Second pattern as aligned by the ground truth."  />
<p class="caption">
Second pattern as aligned by the ground truth.
</p>

While this worked I suppose it is fair to say that our initial pattern
is hardly recognizable. Since we expected this, we planned for a third
kind of pattern in section , that is purely evidence-based. It appears
that, in order to match the ground truths we have at our disposal,
projects register some kind of weak initial peak for the maintenance
activities, that is followed by a somewhat uneventful second and third
interval. Interestingly, the optimization seemed to have used to mostly
straight lines in the Long Stretch phase to model linear declines and
increases. The new Aftermath phase is the longest, so it is clear that
the original pattern and its subdivision into phases is not a good
mapping any longer. Instead of a sharp decline in the Aftermath, we now
see an increase of all variables, without the chance of any decline
before the last observed commit. We will check how this adapted pattern
fares in section .

Pattern III: Averaging the ground truth
---------------------------------------

We can produce a pattern by computing a weighted average over all
available ground truth. As weight, we can use either rater’s score,
their mean or consensus (default).

``` r
gt_weighted_avg <- function(vartype, wtype = c("consensus", "rater.a", "rater.b", "rater.mean")) {
  wtype <- match.arg(wtype)
  gt <- ground_truth[ground_truth[[wtype]] > 0, ]
  wTotal <- sum(gt[[wtype]])
  proj <- gt$project
  weights <- `names<-`(gt[[wtype]], gt$project)
  
  funcs <- lapply(
    project_signals, function(ps) ps[[vartype]]$get0Function())
  
  Vectorize(function(x) {
    val <- 0
    for (p in proj) {
      val <- val + weights[[p]] * funcs[[p]](x)
    }
    val / wTotal
  })
}
```

Now we can easily call above function to produce a weighted average of
each signal:

``` r
p3_avg_signals <- list()

for (vartype in names(weight_vartype)) {
  p3_avg_signals[[vartype]] <- Signal$new(
    name = paste0("p3_avg_", vartype), support = c(0, 1), isWp = TRUE,
    func = gt_weighted_avg(vartype = vartype))
}
```

The 2nd pattern, as derived from the ground truth, is shown in figure .

<img src="fire-drill-technical-report_files/figure-markdown_github/p3-avg-signals-1.png" alt="The third kind of pattern as weighted average over all ground truth."  />
<p class="caption">
The third kind of pattern as weighted average over all ground truth.
</p>

Pattern III (b): Evidence-based
-------------------------------

A third kind of pattern is produced by starting with an empty warping
pattern and having all available ground truth adapt to it. Empty means
that we will start with a flat line located at 0.5 for each variable.
Finally, the parameters are inversed. While we could do this the other
way round, we have two reasons to do it this way, which is the same as
we used for pattern II. First of all if the warping candidate was a
perfectly flat line, it would be very difficult for the gradient to
converge towards some alignment. Secondly, we want to use
equidistantly-spaced boundaries (resulting in equal-length intervals)
and using this approach, we can guarantee the interval lengths. To find
the optimum amount of intervals, we try all values in a certain range
and compute a fit, and then use an information criterion to decide which
of the produced patterns provides the best trade-off between number of
parameters and goodness-of-fit.

The process is the same as for pattern II: Using an instance of
`srBTAW_MultiVartype` that holds one instance of an `srBTAW` per
variable-type. We will choose equidistantly-spaced boundaries over the
WP, and start with just 1 interval, going up to some two-digit number.
The best amount of parameters (intervals) is then determined using the
Akaike Information Criterion (Akaike 1981), which is directly
implemented in `srBTAW`. We either have to use continuous losses or make
sure to **always** use the exact same amount of samples total. The
amount per interval is determined by dividing by the number of
intervals. This is important, as otherwise the information criterion
will not work. We will do a single RSS-loss that covers all intervals.
We will also use an instance of `TimeWarpRegularization` with the
`exint2`-regularizer, as it scales with arbitrary many intervals
(important!). I do not suppose that regularization for the y-values is
needed, so we will not have this. This means that the resulting
objective has just two losses.

For a set of equal-length number of intervals, we will fit such a
multiple variable-type model. This also means we can do this in
parallel. However, models with more intervals and hence more parameters
will considerable take longer during gradient iterations. The more
parameters, the fewer of these models should be fit simultaneously. We
have access to 128-thread machine (of which about 125 thread can be
used). Gradients are computed in parallel as well.

### Preparation

We define a single function that encapsulates the multiple variable-type
model, losses and objectives and returns them, so that we can just fit
them in a loop. The only configurable parameters is the amount of
intervals.

``` r
p3_prepare_mvtypemodel <- function(numIntervals) {
  eps <- sqrt(.Machine$double.eps)
  p3_smv <- srBTAW_MultiVartype$new()
  
  p3_vars <- c("A", "CP", "FREQ", "SCD")
  p3_inst <- list()
  
  # The objective:
  p3_lls <- srBTAW_LossLinearScalarizer$new(
    returnRaw = FALSE, computeParallel = TRUE, gradientParallel = TRUE)
  
  for (name in p3_vars) {
    p3_inst[[name]] <- srBTAW$new(
      # Always includes 0,1 - just as we need it! Works for values >= 1
      theta_b = seq(from = 0, to = 1, by = 1 / numIntervals),
      gamma_bed = c(0, 1, eps),
      lambda = rep(eps, numIntervals),
      begin = 0, end = 1, openBegin = FALSE, openEnd = FALSE,
      useAmplitudeWarping = TRUE,
      # We allow these to be larger; however, the final result should be within [0,1]
      lambda_ymin = rep(-10, numIntervals),
      lambda_ymax = rep( 10, numIntervals),
      isObjectiveLogarithmic = TRUE,
      paramNames = c("v",
        paste0("vtl_", seq_len(length.out = length(p3_vars))),
        paste0("vty_", seq_len(length.out = length(p3_vars)))))
    
    # The WP is a flat line located at 0.5:
    p3_inst[[name]]$setSignal(signal = Signal$new(
      func = function(x) .5, isWp = TRUE, support = c(0, 1), name = paste0("p3_", name)))
    
    # Set the common objective:
    p3_inst[[name]]$setObjective(obj = p3_lls)
    
    # .. and also all the projects' signals:
    for (project in ground_truth[ground_truth$consensus > 0, ]$project) {
      p3_inst[[name]]$setSignal(signal = project_signals[[project]][[name]])
    }
    
    p3_smv$setSrbtaw(varName = name, srbtaw = p3_inst[[name]])
  }

  # We call this there so there are parameters present.
  set.seed(1337 * numIntervals)
  p3_smv$setParams(params =
    `names<-`(x = runif(n = p3_smv$getNumParams()), value = p3_smv$getParamNames()))
  
  for (name in p3_vars) {
    # Add RSS-loss per variable-pair:
    for (project in ground_truth[ground_truth$consensus > 0, ]$project) {
      # The RSS-loss:
      lossRss <- srBTAW_Loss_Rss$new(
        wpName = paste0("p3_", name), wcName = paste(project, name, sep = "_"),
        weight = 1, intervals = seq_len(length.out = numIntervals), continuous = FALSE,
        numSamples = rep(round(5000 / numIntervals), numIntervals), returnRaw = TRUE)
      p3_inst[[name]]$addLoss(loss = lossRss)
      p3_lls$setObjective(
        name = paste(project, name, "rss", sep = "_"), obj = lossRss)
    }
  }
  
  # This has a much higher weight than we had for pattern II
  # because we are using many more samples in the RSS-loss.
  p3_lls$setObjective(name = "p3_reg_exint2", obj = TimeWarpRegularization$new(
    weight = p3_lls$getNumObjectives(), use = "exint2", returnRaw = TRUE,
    wpName = "p3_A", wcName = project_signals$project_1$A$getName(),
    intervals = seq_len(numIntervals)
  )$setSrBtaw(srbtaw = p3_inst$A))
  
  list(smv = p3_smv, lls = p3_lls)
}
```

Now we can compute these in parallel:

``` r
for (numIntervals in c(1:16)) {
  loadResultsOrCompute(
    file = paste0("../results/p3-compute/i_", numIntervals, ".rds"),
    computeExpr =
  {
    p3_vars <- c("A", "CP", "FREQ", "SCD")
    temp <- p3_prepare_mvtypemodel(numIntervals = numIntervals)
    tempf <- tempfile()
    saveRDS(object = temp, file = tempf)
    
    # It does not scale well beyond that.
    cl <- parallel::makePSOCKcluster(min(32, parallel::detectCores()))
    parallel::clusterExport(cl = cl, varlist = list("tempf"))
    
    optR <- doWithParallelClusterExplicit(cl = cl, expr = {
      optimParallel::optimParallel(
        par = temp$smv$getParams(),
        method = "L-BFGS-B",
        lower = c(
          rep(-.Machine$double.xmax, length(p3_vars)), # v_[vartype]
          rep(sqrt(.Machine$double.eps), length(p3_vars)), # vtl
          rep(-.Machine$double.xmax, length(p3_vars) * length(weight_vartype))), # vty for each
        upper = c(
          rep(.Machine$double.xmax, length(p3_vars)),
          rep(1, length(p3_vars)),
          rep(.Machine$double.xmax, length(p3_vars) * length(weight_vartype))),
        fn = function(x) {
          temp <- readRDS(file = tempf)
          temp$smv$setParams(params = x)
          temp$lls$compute0()
        },
        parallel = list(cl = cl, forward = FALSE, loginfo = TRUE)
      )
    })
    list(optR = optR, smv = temp$smv, lls = temp$lls)
  })
}
```

### Finding the best fit

We will load all the results previously computed and compute an
information criterion to compare fits, and then choose the best model.

``` r
p3_params <- NULL
for (tempPath in gtools::mixedsort(
  Sys.glob(paths = paste0(getwd(), "/../results/p3-compute/i*.rds")))
) {
  temp <- readRDS(file = tempPath)
  p3_params <- rbind(p3_params, data.frame(
    numInt = (temp$lls$getNumParams() - 1) / 2,
    numPar = temp$smv$getNumParams(),
    numParSrBTAW = temp$lls$getNumParams(),
    # This AIC would the original one!
    AIC = 2 * temp$smv$getNumParams() - 2 * log(1 / exp(temp$optR$value)),
    # This AIC is based on the number of intervals, not parameters!
    AIC1 = (temp$lls$getNumParams() - 1) - 2 * log(1 / exp(temp$optR$value)),
    # This AIC is based on the amount of parameters per srBTAW instance:
    AIC2 = 2 * temp$lls$getNumParams() - 2 * log(1 / exp(temp$optR$value)),
    logLoss = temp$optR$value,
    loss = exp(temp$optR$value)
  ))
}
```

In table , we show computed fits for various models, where the only
difference is the number of intervals. Each interval comes with two
degrees of freedom: its length and terminal y-translation. Recall that
each computed fit concerns four variables. For example, the first model
with just one interval per variable has nine parameters: All of the
variables share the interval’s length, the first parameter. Then, each
variable has one *v*-parameter, the global y-translation. For each
interval, we have one terminal y-translation. For example, the model
with 7 intervals has 7 + 4 + (4 × 7) = 39 parameters.

We compute the AIC for each fit, which is formulated as in the
following. The parameter *k* is the number of parameters in the model,
i.e., as described, it refers to all the parameters in the
`srBTAW_MultiVartype`-model. The second AIC-alternative uses the
parameter *p* instead, which refers to the number of variables per
`srBTAW`-instance.

$$
\\begin{aligned}
  \\operatorname{AIC}=&\\;2\\times k - 2\\times\\log{(\\mathcal{\\hat{L}})}\\;\\text{, where}\\;\\mathcal{\\hat{L}}\\;\\text{is the maximum log-likelihood of the model,}
  \\\\\[1ex\]
  \\mathcal{\\hat{L}}=&\\;\\frac{1}{\\exp{\\big(\\;\\text{lowest loss of the model}\\;\\big)}}\\;\\text{, since we use logarithmic losses.}
  \\\\\[1em\]
  \\text{The alternatives}&\\;\\operatorname{AIC^1}\\;\\text{and}\\;\\operatorname{AIC^2}\\;\\text{ are defined as:}
  \\\\\[1ex\]
  \\operatorname{AIC^1}=&\\;k-2\\times\\log{(\\mathcal{\\hat{L}})}-1\\;\\text{, which is based on the number of intervals, and}
  \\\\\[1ex\]
  \\operatorname{AIC^2}=&\\;2\\times p - 2\\times\\log{(\\mathcal{\\hat{L}})}\\;\\text{, where}\\;p\\;\\text{is the amount of params per}\\;\\operatorname{srBTAW}\\text{-instance.}
\\end{aligned}
$$

|  numInt|  numPar|  numParSrBTAW|      AIC|    AIC1|    AIC2|  logLoss|       loss|
|-------:|-------:|-------------:|--------:|-------:|-------:|--------:|----------:|
|       1|       9|             3|   34.966|  18.966|  22.966|    8.483|   4831.189|
|       2|      14|             5|   44.052|  20.052|  26.052|    8.026|   3059.055|
|       3|      19|             7|   54.473|  22.473|  30.473|    8.237|   3776.854|
|       4|      24|             9|   63.694|  23.694|  33.694|    7.847|   2557.439|
|       5|      29|            11|   73.985|  25.985|  37.985|    7.992|   2958.219|
|       6|      34|            13|   83.761|  27.761|  41.761|    7.881|   2645.198|
|       7|      39|            15|   94.783|  30.783|  46.783|    8.392|   4410.053|
|       8|      44|            17|  110.695|  38.695|  56.695|   11.347|  84745.003|
|       9|      49|            19|  113.807|  33.807|  53.807|    7.904|   2706.944|
|      10|      54|            21|  123.809|  35.809|  57.809|    7.905|   2709.809|
|      11|      59|            23|  134.167|  38.167|  62.167|    8.084|   3240.653|
|      12|      64|            25|  144.350|  40.350|  66.350|    8.175|   3551.062|
|      13|      69|            27|  153.878|  41.878|  69.878|    7.939|   2804.353|
|      14|      74|            29|  164.122|  44.122|  74.122|    8.061|   3168.120|
|      15|      79|            31|  174.284|  46.284|  78.284|    8.142|   3435.957|
|      16|      84|            33|  184.305|  48.305|  82.305|    8.153|   3472.327|

Comparing the results from table , it appears that no matter how we
define the AIC, it is increasing with the number of parameters, and it
does so faster than the loss reduces. So, picking a model by AIC is not
terribly useful, as the results suggest we would to go with the
1-interval model. The model with the lowest loss is the one with 4
intervals.

### Create pattern from best fit

This is the same process as for pattern II, as the parameters need
inversion. We will reconstruct the warped signals according to the
inversed parameters to produce the third pattern. According to the
overview above, the best model (lowest loss, **not** AIC) is the one
with **4** intervals. Its parameters are the following:

``` r
# Let's first define a function that inverses the params and reconstructs the pattern.
p3_pattern_from_fit <- function(whichNumIntervals) {
  res <- list()
  p3_i <- readRDS(file = paste0(getwd(), "/../results/p3-compute/i_", whichNumIntervals, ".rds"))
  
  # FitResult:
  fr <- FitResult$new("foo")
  fr$fromOptimParallel(optR = p3_i$optR)
  res$fr <- fr
  
  # Inversion:
  lambda <- p3_i$smv$.__enclos_env__$private$instances$A$.__enclos_env__$private$instances$`p3_A|project_1_A`$getLambda()
  p3_i_varthetaL <- p3_i$optR$par[grepl(pattern = "^vtl_", x = names(p3_i$optR$par))]
  for (q in seq_len(length.out = whichNumIntervals)) {
    if (p3_i_varthetaL[q] < lambda[q]) {
      p3_i_varthetaL[q] <- lambda[q]
    }
  }
  p3_i_varthetaL <- p3_i_varthetaL / sum(p3_i_varthetaL)
  
  p3_i_thetaB <- c(0)
  for (idx in seq_len(length.out = length(p3_i_varthetaL))) {
    p3_i_thetaB <- c(p3_i_thetaB, sum(p3_i_varthetaL[1:idx]))
  }
  p3_i_thetaB[length(p3_i_thetaB)] <- 1 # numeric stability
  
  p3_i_varthetaL
  p3_i_thetaB
  res$varthetaL <- p3_i_varthetaL
  res$thetaB <- p3_i_thetaB
  
  # Signals:
  p3_i_numInt <- length(p3_i_varthetaL)
  p3_i_signals <- list()
  
  for (vartype in names(weight_vartype)) {
    emptySig <- Signal$new(
      isWp = TRUE, # does not matter here
        func = function(x) .5, support = c(0, 1), name = paste0("p3_", vartype))
    
    temp <- SRBTWBAW$new(
      theta_b = unname(p3_i_thetaB), gamma_bed = c(0, 1, 0),
      wp = emptySig$get0Function(), wc = emptySig$get0Function(),
      lambda = rep(0, p3_i_numInt), begin = 0, end = 1,
      lambda_ymin = rep(0, p3_i_numInt), lambda_ymax = rep(1, p3_i_numInt))
    
    # Recall that originally we used equidistantly-spaced boundaries:
    temp$setParams(vartheta_l = rep(1 / p3_i_numInt, p3_i_numInt),
                   v = -1 * p3_i$optR$par[paste0("v_", vartype)],
                   vartheta_y = -1 * p3_i$optR$par[paste0("vty_", 1:p3_i_numInt, "_", vartype)])
    
    p3_i_signals[[vartype]] <- Signal$new(
      name = paste0("p3_", vartype), support = c(0, 1), isWp = TRUE, func = Vectorize(temp$M))
  }
  res$signals <- p3_i_signals
  
  # Data:
  temp <- NULL
  for (vartype in names(weight_vartype)) {
    f <- p3_i_signals[[vartype]]$get0Function()
    x <- seq(from = 0, to = 1, length.out = 1e3)
    y <- f(x)
    
    temp <- rbind(temp, data.frame(
      x = x,
      y = y,
      t = vartype,
      numInt = whichNumIntervals
    ))
  }
  res$data <- temp
  res
}
```

``` r
p3_best <- readRDS(file = paste0(getwd(), "/../results/p3-compute/i_",
                     p3_params[which.min(p3_params$loss), ]$numInt, ".rds"))
p3_best$optR$par
```

    ##           v_A          v_CP        v_FREQ         v_SCD         vtl_1 
    ##  0.4753518248  0.4449292077  0.3969243406 -0.3425136575  0.0486178486 
    ##         vtl_2         vtl_3         vtl_4       vty_1_A      vty_1_CP 
    ##  0.0160263762  0.3063930774  0.9919315307 -0.0085322067  0.0130265463 
    ##    vty_1_FREQ     vty_1_SCD       vty_2_A      vty_2_CP    vty_2_FREQ 
    ##  0.0423862173  0.0196908380 -0.0055657708 -0.0293088290 -0.0404751674 
    ##     vty_2_SCD       vty_3_A      vty_3_CP    vty_3_FREQ     vty_3_SCD 
    ##  0.0401942364 -0.0898785537 -0.0488678587 -0.1533263719  0.0002213948 
    ##       vty_4_A      vty_4_CP    vty_4_FREQ     vty_4_SCD 
    ##  0.0673046097 -0.1805227900 -0.1085880858 -0.0666110009

First we have to inverse the parameters before we can reconstruct the
signals:

``` r
p3_best_varthetaL <- p3_best$optR$par[
  grepl(pattern = "^vtl_", x = names(p3_best$optR$par))]
p3_best_varthetaL <- p3_best_varthetaL / sum(p3_best_varthetaL)

p3_best_thetaB <- c(0)
for (idx in seq_len(length.out = length(p3_best_varthetaL))) {
  p3_best_thetaB <- c(p3_best_thetaB, sum(p3_best_varthetaL[1:idx]))
}
p3_best_thetaB[length(p3_best_thetaB)] <- 1 # numeric stability

p3_best_varthetaL
```

    ##      vtl_1      vtl_2      vtl_3      vtl_4 
    ## 0.03567055 0.01175843 0.22479830 0.72777272

``` r
p3_best_thetaB
```

    ## [1] 0.00000000 0.03567055 0.04742898 0.27222728 1.00000000

``` r
p3_best_numInt <- length(p3_best_varthetaL)
p3_signals <- list()

for (vartype in names(weight_vartype)) {
  emptySig <- Signal$new(
    isWp = TRUE, # does not matter here
      func = function(x) .5, support = c(0, 1), name = paste0("p3_", vartype))
  
  temp <- SRBTWBAW$new(
    theta_b = unname(p3_best_thetaB), gamma_bed = c(0, 1, 0),
    wp = emptySig$get0Function(), wc = emptySig$get0Function(),
    lambda = rep(0, p3_best_numInt), begin = 0, end = 1,
    lambda_ymin = rep(0, p3_best_numInt), lambda_ymax = rep(1, p3_best_numInt))
  
  # Recall that originally we used equidistantly-spaced boundaries:
  temp$setParams(vartheta_l = rep(1 / p3_best_numInt, p3_best_numInt),
                 v = -1 * p3_best$optR$par[paste0("v_", vartype)],
                 vartheta_y = -1 * p3_best$optR$par[paste0("vty_", 1:p3_best_numInt, "_", vartype)])
  
  p3_signals[[vartype]] <- Signal$new(
    name = paste0("p3_", vartype), support = c(0, 1), isWp = TRUE, func = Vectorize(temp$M))
}
```

The 2nd pattern, as derived from the ground truth, is shown in figure .

<img src="fire-drill-technical-report_files/figure-markdown_github/p3-signals-1.png" alt="Pattern type III (b) pattern as aligned by the ground truth only."  />
<p class="caption">
Pattern type III (b) pattern as aligned by the ground truth only.
</p>

Let’s show all computed patterns in a grid:

In figure we can clearly observe how the pattern evolves with growing
number or parameters. Almost all patterns with sufficiently many degrees
of freedom have some crack at about one quarter of the projects’ time, a
second crack is observed at about three quarter’s time. In all patterns,
it appears that adaptive activities are the least common. All patterns
started with randomized coefficients, and something must have gone wrong
for pattern 8. From five and more intervals we can observe growing
similarities with the weighted-average pattern, although it never comes
really close. Even though we used a timewarp-regularizer with high
weight, we frequently get extreme intervals.

<img src="fire-drill-technical-report_files/figure-markdown_github/p3-all-1.png" alt="Computed pattern by number of intervals."  />
<p class="caption">
Computed pattern by number of intervals.
</p>

In figure we can clearly see that all but the eighth pattern converged
nicely (this was already visible in ). The loss is logarithmic, so the
progress is rather substantial. For example, going from
log (14) ≈ 1.2*e*6 to log (8) ≈ 3*e*3 is a reduction by 3 (!) orders of
magnitude.

<img src="fire-drill-technical-report_files/figure-markdown_github/p3b-all-fr-1.png" alt="Losses for all computed pattern by number of intervals."  />
<p class="caption">
Losses for all computed pattern by number of intervals.
</p>

Scoring of projects
===================

The true main-purpose of our work is to take a pattern and check it
against any project, with the goal of obtaining a score, or
goodness-of-match so that we can determine if the AP in the pattern is
present in the project. In the previous sections we have introduced a
number of patterns that we are going to apply here.

How it works: Given some pattern that consists of one or arbitrary many
signals, the pattern is added to a single instance of `srBTAW` as
**Warping Pattern**. The project’s signals are added as **Warping
Candidates** to the same instance.

To compute a score, we need to define how to measure the distance
between the WP and the WC (between each pair of signals and each
interval). In the notebooks for sr-BTW we have previously defined some
suitable losses with either **global** or **local** finite upper bounds.
Currently, the Jensen–Shannon divergence (JSD), as well as the
ratio-metrics (correlation, arc-lengths) have global upper bounds. For
the JSD, it is ln 2. Losses with local finite upper bound are, for
example, the area between curves, the residual sum of squares, the
Euclidean distance etc., basically any metric that has a limit within
the rectangle demarcated by one or more intervals. For some of the
patterns, we have used a combination of such losses with local bounds.
In general, it is not necessary to fit a pattern with the same kinds of
losses that are later on used for scoring, but it is recommended to
avoid confusing may results.

Scoring mechanisms
------------------

For scoring a single project, we first warp it to the pattern, then we
measure the remaining distance. We only do time-warping of the projects
to the pattern. We could compute a score for each interval. However, the
ground truth does not yield this, so we will only compute a scores for
entire signals, i.e., over all intervals. Once aligned, computing scores
is cheap, so we will try a variety of scores and see what works best.

``` r
score_variable_alignment <- function(
  alignment, patternName, projectName,
  use = c("area", "corr", "jsd", "kl", "arclen", "sd", "var", "mae",
          "rmse", "RMS", "Kurtosis", "Peak", "ImpulseFactor")
) {
  use <- match.arg(use)
  sigs <- alignment$getAllSignals()
  vartypes <- names(weight_vartype)
  
  scores <- c()
  for (vartype in vartypes) {
    inst <- alignment$getInstance(
      wpName = paste0(patternName, "_", vartype),
      wcName = paste0(projectName, "_", vartype))
    supp <- c(inst$getTb_q(q = 1), inst$getTe_q(q = max(inst$getQ())))
    stopifnot(supp[1] == 0 && supp[2] == 1)
    # Now the two functions/signals: The first is the unchanged WP,
    # the second is the warped M function of the SRBTW/SRBTWBAW!
    f1 <- inst$getWP()
    f2 <- Vectorize(inst$M)
    temp <- NA_real_
    
    if (use == "area") {
      temp <- area_diff_2_functions_score()(f1 = f1, f2 = f2)
    } else if (use == "corr") {
      tempsd <- stat_diff_2_functions(f1 = f1, f2 = f2)
      temp <- stats::cor(x = tempsd$dataF1, y = tempsd$dataF2)
      if (is.na(temp)) {
        temp <- 0
      }
      temp <- (temp + 1) / 2 # move from [-1,1] \to [0,1]
    } else if (use == "jsd") {
      temp <- stat_diff_2_functions_symmetric_JSD_score()(f1 = f1, f2 = f2)
    } else if (use == "kl") {
      temp <- stat_diff_2_functions_symmetric_KL_sampled(f1 = f1, f2 = f2)
      temp <- stat_diff_2_functions_mutual_information_score(
        sensitivityExponent = 5)(f1 = f1, f2 = f2)$value
    } else if (use == "arclen") {
      temp <- stat_diff_2_functions_arclen_score()(f1 = f1, f2 = f2)
    } else if (use %in% c("RMS", "Kurtosis", "Peak", "ImpulseFactor")) {
      temp <- stat_diff_2_functions_signals_score(use = use)(f1 = f1, f2 = f2)
    } else if (use %in% c("sd", "var", "mae", "rmse")) {
      temp <- stat_diff_2_functions_sd_var_mae_rmse_score(use = use)(f1 = f1, f2 = f2)
    } else stop(paste0("Don't know score ", use))
    
    scores <- c(scores, `names<-`(c(temp), vartype))
  }
  scores
}
```

We define a parallelized function to compute all scores of a project:

``` r
compute_all_scores <- function(alignment, patternName) {
  useScores <- c("area", "corr", "jsd", "kl", "arclen", "sd", "var",
                 "mae", "rmse", "RMS", "Kurtosis", "Peak", "ImpulseFactor")
  
  #`rownames<-`(
  doWithParallelCluster(numCores = length(alignment), expr = {
    foreach::foreach(
      projectName = names(alignment),
      .inorder = TRUE,
      .combine = rbind,
      .export = c("score_variable_alignment", "weight_vartype")
    ) %dopar% {
      source("./common-funcs.R")
      source("../models/modelsR6.R")
      source("../models/SRBTW-R6.R")
      
      scores <- c()
      for (score in useScores) {
        temp <- score_variable_alignment(
          patternName = patternName, projectName = projectName,
          alignment = alignment[[projectName]], use = score)
        scores <- c(scores, `names<-`(
          c(mean(temp), prod(temp)), c(paste0(score, c("_m", "_p")))))
      }
      `colnames<-`(matrix(data = scores, nrow = 1), names(scores))
    }
  })#, sort(names(alignment)))
}
```

We also need to define a function for warping a project to the pattern:

``` r
#' This function creates a new instance of \code{srBTAW} and adds all
#' the signals of pattern 1 and the given project to it. Then it creates
#' a single RSS-loss per variable pair and adds it to the model and a
#' linear scalarizer that is used as objective.
time_warp_project <- function(pattern, project, thetaB = c(0, fd_data_boundaries, 1)) {
  stopifnot(is.logical(all.equal(names(pattern), names(project))))
  numInt <- length(thetaB) - 1
  
  obj <- srBTAW_LossLinearScalarizer$new(
    computeParallel = TRUE, gradientParallel = TRUE, returnRaw = FALSE)
  
  inst <- srBTAW$new(
    theta_b = thetaB,
    gamma_bed = c(0, 1, sqrt(.Machine$double.eps)),
    lambda = rep(sqrt(.Machine$double.eps), numInt),
    begin = 0, end = 1, openBegin = FALSE, openEnd = FALSE,
    useAmplitudeWarping = FALSE)
  inst$setParams(params = `names<-`(rep(1/numInt, numInt), inst$getParamNames()))
  
  for (vartype in names(pattern)) {
    inst$setSignal(signal = pattern[[vartype]])
    inst$setSignal(signal = project[[vartype]])
    
    loss <- srBTAW_Loss_Rss$new(
      intervals = seq_len(length.out = numInt), returnRaw = TRUE, # NOT logarithmic!
      wpName = pattern[[vartype]]$getName(), wcName = project[[vartype]]$getName(),
      weight = 1, continuous = FALSE, numSamples = rep(1e3, numInt))
    inst$addLoss(loss = loss)
    obj$setObjective(name = paste("rss", vartype, sep = "_"), obj = loss)
  }
  
  # Add a time-warp regularizer:
  reg <- TimeWarpRegularization$new(
    weight = 1/2, use = "exint2", wpName = pattern$A$getName(), wcName = project$A$getName(),
    returnRaw = TRUE, intervals = seq_len(length.out = length(names(project))))
  inst$addLoss(loss = reg)
  obj$setObjective(name = "reg_exint2", obj = reg)

  inst$setObjective(obj = obj)
}
```

Pattern I
---------

First we compute the alignment for all projects, then all scores.

``` r
library(foreach)

p1_align <- loadResultsOrCompute(file = "../results/p1_align.rds", computeExpr = {
  # Let's compute all projects in parallel!
  cl <- parallel::makePSOCKcluster(length(project_signals))
  unlist(doWithParallelClusterExplicit(cl = cl, expr = {
    foreach::foreach(
      projectName = names(project_signals),
      .inorder = FALSE,
      .packages = c("parallel")
    ) %dopar% {
      source("./common-funcs.R")
      source("../models/modelsR6.R")
      source("../models/SRBTW-R6.R")
      
      # There are 5 objectives that can be computed in parallel!
      cl_nested <- parallel::makePSOCKcluster(5)
      `names<-`(list(doWithParallelClusterExplicit(cl = cl_nested, expr = {
        temp <- time_warp_project(
          pattern = p1_signals, project = project_signals[[projectName]])
        temp$fit(verbose = TRUE)
        temp # return the instance, it includes the FitResult
      })), projectName)
    }
  }))
})
```

``` r
p1_scores <- loadResultsOrCompute(file = "../results/p1_scores.rds", computeExpr = {
  as.data.frame(compute_all_scores(alignment = p1_align, patternName = "p1"))
})
```

Recall that we are obtaining scores for each interval. To aggregate them
we build the product and the mean in the following table, there is no
weighing applied.

|                  |   pr\_1|      pr\_2|  pr\_3|  pr\_4|   pr\_5|     pr\_6|     pr\_7|  pr\_8|  pr\_9|
|:-----------------|-------:|----------:|------:|------:|-------:|---------:|---------:|------:|------:|
| area\_m          |    0.80|       0.81|   0.89|   0.89|    0.83|      0.83|      0.77|   0.84|   0.82|
| area\_p          |    0.40|       0.42|   0.63|   0.63|    0.48|      0.48|      0.34|   0.50|   0.45|
| corr\_m          |    0.50|       0.64|   0.67|   0.58|    0.59|      0.52|      0.55|   0.54|   0.54|
| corr\_p          |    0.05|       0.15|   0.19|   0.11|    0.10|      0.03|      0.07|   0.08|   0.06|
| jsd\_m           |    0.28|       0.34|   0.46|   0.43|    0.40|      0.37|      0.29|   0.34|   0.32|
| jsd\_p           |    0.00|       0.01|   0.03|   0.02|    0.02|      0.01|      0.00|   0.01|   0.01|
| kl\_m            |    5.74|      85.40|   4.19|   4.35|    6.22|    233.32|     37.26|   3.32|   3.95|
| kl\_p            |  107.94|  169731.21|  62.15|  74.76|  203.24|  35022.00|  10435.03|  38.18|  14.37|
| arclen\_m        |    0.47|       0.52|   0.43|   0.51|    0.57|      0.85|      0.54|   0.53|   0.50|
| arclen\_p        |    0.04|       0.07|   0.03|   0.06|    0.10|      0.50|      0.06|   0.08|   0.06|
| sd\_m            |    0.65|       0.70|   0.77|   0.74|    0.76|      0.75|      0.69|   0.70|   0.68|
| sd\_p            |    0.18|       0.24|   0.34|   0.30|    0.32|      0.30|      0.21|   0.23|   0.21|
| var\_m           |    0.87|       0.91|   0.94|   0.93|    0.94|      0.93|      0.89|   0.90|   0.89|
| var\_p           |    0.58|       0.67|   0.79|   0.75|    0.77|      0.75|      0.63|   0.66|   0.64|
| mae\_m           |    0.80|       0.81|   0.89|   0.89|    0.83|      0.83|      0.77|   0.84|   0.82|
| mae\_p           |    0.40|       0.42|   0.62|   0.63|    0.48|      0.48|      0.34|   0.50|   0.45|
| rmse\_m          |    0.74|       0.76|   0.83|   0.81|    0.78|      0.80|      0.71|   0.76|   0.76|
| rmse\_p          |    0.30|       0.33|   0.48|   0.43|    0.37|      0.40|      0.24|   0.33|   0.32|
| RMS\_m           |    0.66|       0.73|   0.50|   0.60|    0.52|      0.76|      0.58|   0.69|   0.50|
| RMS\_p           |    0.14|       0.24|   0.03|   0.06|    0.04|      0.28|      0.06|   0.18|   0.04|
| Kurtosis\_m      |    0.40|       0.49|   0.39|   0.31|    0.24|      0.39|      0.35|   0.35|   0.13|
| Kurtosis\_p      |    0.00|       0.00|   0.00|   0.00|    0.00|      0.00|      0.00|   0.00|   0.00|
| Peak\_m          |    0.62|       0.62|   0.53|   0.63|    0.65|      0.61|      0.61|   0.63|   0.44|
| Peak\_p          |    0.09|       0.10|   0.04|   0.06|    0.12|      0.09|      0.06|   0.09|   0.03|
| ImpulseFactor\_m |    0.60|       0.59|   0.54|   0.66|    0.48|      0.61|      0.57|   0.54|   0.64|
| ImpulseFactor\_p |    0.11|       0.11|   0.05|   0.11|    0.05|      0.08|      0.06|   0.06|   0.16|

``` r
corr <- stats::cor(ground_truth$consensus, p1_scores)[1, ]

if (interactive()) {
  corr
} else {
  perCol <- ceiling(length(corr) / 3)
  temp <- data.frame(matrix(ncol = 6, nrow = perCol))
  for (idx in 1:3) {
    off <- (idx - 1) * perCol
    temp[, idx * 2 - 1] <- names(corr)[(1 + off):(perCol + off)]
    temp[, idx * 2] <- corr[(1 + off):(perCol + off)]
  }
  colnames(temp) <- rep(c("Score", "Value"), 2)
  knitr::kable(
    x = temp,
    booktabs = TRUE,
    caption = "Correlation of the ground truth with all other scores for pattern I.",
    label = "p1-corr")
}
```

| Score     |       Value| Score     |       Value| NA               |          NA|
|:----------|-----------:|:----------|-----------:|:-----------------|-----------:|
| area\_m   |   0.6115889| arclen\_p |  -0.1684936| RMS\_m           |  -0.5768954|
| area\_p   |   0.6436132| sd\_m     |   0.3796681| RMS\_p           |  -0.6222672|
| corr\_m   |   0.2684177| sd\_p     |   0.3934965| Kurtosis\_m      |  -0.3730592|
| corr\_p   |   0.2605112| var\_m    |   0.3842355| Kurtosis\_p      |  -0.6263601|
| jsd\_m    |   0.5388574| var\_p    |   0.3900069| Peak\_m          |  -0.4524193|
| jsd\_p    |   0.5938024| mae\_m    |   0.6101791| Peak\_p          |  -0.7678158|
| kl\_m     |  -0.2513813| mae\_p    |   0.6423594| ImpulseFactor\_m |   0.4717221|
| kl\_p     |  -0.4087868| rmse\_m   |   0.4837649| ImpulseFactor\_p |   0.2300525|
| arclen\_m |  -0.2546068| rmse\_p   |   0.5246423| NA               |          NA|

Let’s show a correlation matrix in figure :

<img src="fire-drill-technical-report_files/figure-markdown_github/p1-corr-mat-1.png" alt="Correlation matrix for scores using pattern I."  />
<p class="caption">
Correlation matrix for scores using pattern I.
</p>

We appear to have mostly strongly negative correlations – note that the
measure using Kullback-Leibler is a divergence, not a similarity. Area
and RMSE have a strong negative correlation, which suggests that
whenever their score is high, the ground truth’s consensus score is low.
A high score for area or RMSE however means, that the distance between
the signals is comparatively low, so we should have a good alignment, so
how is this explained then?

Going back to table , we will notice that projects 2, 3 and 5 have
somewhat similar courses for their variables, yet their consensus scores
are 0, 6 and 1, respectively. In other words, only project 3 has had a
Fire Drill. We can hence conclude that the visual distance and the
ground truth consensus are **not** proportional (at least not for the
variables we chose to model). The visual similarity of a project to a
pattern is just that; the score quantifies the deviation from the
pattern, but it does not necessarily correlate with the ground truth.
This however was our underlying assumption all along, hence the initial
pattern. We deliberately chose to design patterns *without*
investigating any of the projects. Also, while we had access to the
projects for some time now, the ground truth became available only very
recently, after all modeling was done.

Nevertheless, it does not hurt to check out patterns II and III, as we
would like to achieve better matches. Eventually, the goal with this
approach is to improve correlations and to get more accurate scores. The
final stage then could be to compute, for example, a weighted consensus,
based on the projects that we have, or to create a linear model that can
regress to a value close to the ground truth by considering all the
different scores.

Pattern II
----------

The second pattern was produced by having it warp to all the ground
truths simultaneously, using their weight.

``` r
library(foreach)

p2_align <- loadResultsOrCompute(file = "../results/p2_align.rds", computeExpr = {
  # Let's compute all projects in parallel!
  cl <- parallel::makePSOCKcluster(length(project_signals))
  unlist(doWithParallelClusterExplicit(cl = cl, expr = {
    foreach::foreach(
      projectName = names(project_signals),
      .inorder = FALSE,
      .packages = c("parallel")
    ) %dopar% {
      source("./common-funcs.R")
      source("../models/modelsR6.R")
      source("../models/SRBTW-R6.R")
      
      # There are 5 objectives that can be computed in parallel!
      cl_nested <- parallel::makePSOCKcluster(5)
      `names<-`(list(doWithParallelClusterExplicit(cl = cl_nested, expr = {
        temp <- time_warp_project(
          pattern = p2_signals, project = project_signals[[projectName]])
        temp$fit(verbose = TRUE)
        temp # return the instance, it includes the FitResult
      })), projectName)
    }
  }))
})
```

``` r
p2_scores <- loadResultsOrCompute(file = "../results/p2_scores.rds", computeExpr = {
  #as.data.frame(
  compute_all_scores(alignment = p2_align, patternName = "p2")
  #)
})
```

|                  |      pr\_1|       pr\_2|    pr\_3|         pr\_4|      pr\_5|  pr\_6|   pr\_7|      pr\_8|     pr\_9|
|:-----------------|----------:|-----------:|--------:|-------------:|----------:|------:|-------:|----------:|---------:|
| area\_m          |       0.92|        0.92|     0.90|  9.000000e-01|       0.89|   0.87|    0.86|       0.90|      0.91|
| area\_p          |       0.72|        0.71|     0.65|  6.400000e-01|       0.61|   0.56|    0.55|       0.65|      0.68|
| corr\_m          |       0.69|        0.61|     0.55|  6.800000e-01|       0.58|   0.52|    0.59|       0.69|      0.58|
| corr\_p          |       0.22|        0.13|     0.09|  2.100000e-01|       0.10|   0.05|    0.12|       0.21|      0.11|
| jsd\_m           |       0.48|        0.42|     0.36|  4.100000e-01|       0.35|   0.29|    0.34|       0.34|      0.38|
| jsd\_p           |       0.05|        0.03|     0.02|  2.000000e-02|       0.01|   0.00|    0.01|       0.01|      0.02|
| kl\_m            |      47.29|      126.99|     9.78|  1.064208e+05|      55.98|  20.30|    4.75|      91.93|     52.51|
| kl\_p            |  152516.09|  4876412.13|  1607.26|  5.664986e+10|  425504.55|   9.55|  188.76|  185115.00|  30601.09|
| arclen\_m        |       0.53|        0.54|     0.48|  5.700000e-01|       0.59|   0.73|    0.51|       0.58|      0.57|
| arclen\_p        |       0.05|        0.07|     0.04|  6.000000e-02|       0.11|   0.24|    0.04|       0.09|      0.08|
| sd\_m            |       0.86|        0.83|     0.81|  8.400000e-01|       0.81|   0.76|    0.76|       0.80|      0.82|
| sd\_p            |       0.54|        0.48|     0.44|  5.000000e-01|       0.43|   0.31|    0.32|       0.42|      0.45|
| var\_m           |       0.98|        0.97|     0.96|  9.700000e-01|       0.96|   0.93|    0.94|       0.96|      0.97|
| var\_p           |       0.92|        0.89|     0.86|  9.000000e-01|       0.85|   0.75|    0.76|       0.85|      0.87|
| mae\_m           |       0.92|        0.92|     0.90|  9.000000e-01|       0.89|   0.87|    0.86|       0.90|      0.91|
| mae\_p           |       0.72|        0.71|     0.64|  6.400000e-01|       0.61|   0.56|    0.55|       0.65|      0.68|
| rmse\_m          |       0.89|        0.88|     0.86|  8.600000e-01|       0.85|   0.81|    0.81|       0.85|      0.87|
| rmse\_p          |       0.62|        0.59|     0.55|  5.600000e-01|       0.53|   0.42|    0.42|       0.53|      0.57|
| RMS\_m           |       0.62|        0.66|     0.71|  6.500000e-01|       0.83|   0.62|    0.58|       0.64|      0.83|
| RMS\_p           |       0.11|        0.14|     0.21|  1.500000e-01|       0.45|   0.14|    0.11|       0.13|      0.46|
| Kurtosis\_m      |       0.42|        0.29|     0.22|  3.000000e-01|       0.33|   0.11|    0.30|       0.39|      0.37|
| Kurtosis\_p      |       0.01|        0.00|     0.00|  0.000000e+00|       0.01|   0.00|    0.00|       0.00|      0.00|
| Peak\_m          |       0.66|        0.52|     0.54|  5.300000e-01|       0.70|   0.53|    0.64|       0.66|      0.60|
| Peak\_p          |       0.15|        0.06|     0.07|  6.000000e-02|       0.23|   0.07|    0.14|       0.14|      0.11|
| ImpulseFactor\_m |       0.66|        0.63|     0.57|  5.100000e-01|       0.62|   0.61|    0.77|       0.55|      0.55|
| ImpulseFactor\_p |       0.17|        0.13|     0.07|  2.000000e-02|       0.09|   0.08|    0.34|       0.06|      0.06|

The correlation of just the ground truth with all scores is in table .

| Score     |       Value| Score     |       Value| Score            |       Value|
|:----------|-----------:|:----------|-----------:|:-----------------|-----------:|
| area\_m   |  -0.0890189| arclen\_p |  -0.2557238| RMS\_m           |   0.1676109|
| area\_p   |  -0.1069064| sd\_m     |   0.1285950| RMS\_p           |   0.1289162|
| corr\_m   |  -0.0764780| sd\_p     |   0.1292517| Kurtosis\_m      |  -0.2078491|
| corr\_p   |  -0.0362952| var\_m    |   0.1332358| Kurtosis\_p      |  -0.3601837|
| jsd\_m    |   0.0077401| var\_p    |   0.1339373| Peak\_m          |  -0.4622136|
| jsd\_p    |  -0.1128088| mae\_m    |  -0.0908586| Peak\_p          |  -0.4417711|
| kl\_m     |   0.6724389| mae\_p    |  -0.1087002| ImpulseFactor\_m |  -0.4271623|
| kl\_p     |   0.6729587| rmse\_m   |   0.0336848| ImpulseFactor\_p |  -0.2815468|
| arclen\_m |  -0.1780968| rmse\_p   |   0.0189570| NA               |          NA|

The correlation matrix looks as in figure .

<img src="fire-drill-technical-report_files/figure-markdown_github/p2-corr-mat-1.png" alt="Correlation matrix for scores using pattern II."  />
<p class="caption">
Correlation matrix for scores using pattern II.
</p>

With the second pattern we get much stronger correlations on average,
meaning that the alignment of each project to the second pattern is
better. While some correlations with the ground truth remain similar, we
get strong positive correlations with the signal-measures
*Impulse-factor* (-0.427) and *Peak* (-0.462) (Xi, Sun, and Krishnappa
2000). This however is explained by the double warping: First the
pattern II was produced by time- and amplitude-warping it to the ground
truth. Then, each project was time-warped to the that pattern. Simply
put, this should result in some good alignment of all of the signals’
peaks, explaining the high correlations.

Pattern II (without alignment)
------------------------------

Since pattern II was computed such that it warps to all projects, it
already should correct for time- and amplitude warping. This gives us
incentive to compute scores against the non-aligned projects:

``` r
p2_no_align <- list()

for (project in ground_truth$project) {
  temp <- p2_align[[project]]$clone()
  temp$setParams(params = `names<-`(p2_ex_varthetaL, temp$getParamNames()))
  p2_no_align[[project]] <- temp
}
```

``` r
p2_no_scores <- loadResultsOrCompute(file = "../results/p2_no_scores.rds", computeExpr = {
  as.data.frame(compute_all_scores(alignment = p2_no_align, patternName = "p2"))
})
```

|                  |    pr\_1|     pr\_2|      pr\_3|          pr\_4|   pr\_5|       pr\_6|     pr\_7|   pr\_8|     pr\_9|
|:-----------------|--------:|---------:|----------:|--------------:|-------:|-----------:|---------:|-------:|---------:|
| area\_m          |     0.78|      0.88|       0.89|           0.88|    0.85|        0.86|      0.89|    0.83|      0.90|
| area\_p          |     0.36|      0.61|       0.63|           0.59|    0.52|        0.53|      0.62|    0.48|      0.64|
| corr\_m          |     0.57|      0.68|       0.60|           0.61|    0.55|        0.44|      0.55|    0.62|      0.62|
| corr\_p          |     0.10|      0.21|       0.10|           0.13|    0.09|        0.03|      0.08|    0.13|      0.14|
| jsd\_m           |     0.32|      0.42|       0.36|           0.33|    0.31|        0.26|      0.31|    0.35|      0.37|
| jsd\_p           |     0.01|      0.02|       0.01|           0.01|    0.01|        0.00|      0.01|    0.01|      0.02|
| kl\_m            |    14.64|     23.24|     112.70|         764.58|    4.29|      597.19|     30.19|    5.43|     21.92|
| kl\_p            |  2924.94|  10622.94|  363270.85|  480249909\.36|  139.32|  2771576.46|  49009.42|  249.27|  13800.41|
| arclen\_m        |     0.76|      0.55|       0.55|           0.58|    0.51|        0.57|      0.47|    0.55|      0.60|
| arclen\_p        |     0.31|      0.07|       0.08|           0.06|    0.04|        0.08|      0.04|    0.08|      0.12|
| sd\_m            |     0.73|      0.79|       0.78|           0.79|    0.73|        0.75|      0.79|    0.78|      0.81|
| sd\_p            |     0.27|      0.37|       0.37|           0.39|    0.28|        0.30|      0.38|    0.35|      0.43|
| var\_m           |     0.92|      0.95|       0.95|           0.95|    0.92|        0.93|      0.95|    0.94|      0.96|
| var\_p           |     0.71|      0.80|       0.81|           0.82|    0.71|        0.74|      0.81|    0.79|      0.86|
| mae\_m           |     0.78|      0.88|       0.89|           0.88|    0.85|        0.86|      0.89|    0.83|      0.90|
| mae\_p           |     0.36|      0.61|       0.63|           0.59|    0.52|        0.53|      0.62|    0.48|      0.64|
| rmse\_m          |     0.72|      0.83|       0.84|           0.84|    0.79|        0.81|      0.85|    0.79|      0.86|
| rmse\_p          |     0.26|      0.47|       0.50|           0.49|    0.38|        0.42|      0.51|    0.38|      0.55|
| RMS\_m           |     0.49|      0.63|       0.59|           0.60|    0.54|        0.63|      0.65|    0.64|      0.75|
| RMS\_p           |     0.06|      0.14|       0.09|           0.10|    0.08|        0.15|      0.16|    0.13|      0.32|
| Kurtosis\_m      |     0.13|      0.35|       0.21|           0.17|    0.23|        0.23|      0.17|    0.39|      0.33|
| Kurtosis\_p      |     0.00|      0.00|       0.00|           0.00|    0.00|        0.00|      0.00|    0.00|      0.00|
| Peak\_m          |     0.53|      0.66|       0.52|           0.53|    0.64|        0.60|      0.54|    0.66|      0.70|
| Peak\_p          |     0.07|      0.15|       0.06|           0.06|    0.14|        0.11|      0.07|    0.14|      0.23|
| ImpulseFactor\_m |     0.64|      0.82|       0.65|           0.49|    0.75|        0.51|      0.57|    0.60|      0.64|
| ImpulseFactor\_p |     0.13|      0.44|       0.16|           0.04|    0.30|        0.04|      0.08|    0.11|      0.14|

The correlation of just the ground truth with all scores is in table .

| Score     |       Value| Score     |       Value| Score            |       Value|
|:----------|-----------:|:----------|-----------:|:-----------------|-----------:|
| area\_m   |   0.5183212| arclen\_p |  -0.1900343| RMS\_m           |   0.2211204|
| area\_p   |   0.5317444| sd\_m     |   0.5396011| RMS\_p           |   0.1957347|
| corr\_m   |   0.0640167| sd\_p     |   0.5611052| Kurtosis\_m      |  -0.4022877|
| corr\_p   |  -0.0668462| var\_m    |   0.5762125| Kurtosis\_p      |  -0.1080604|
| jsd\_m    |  -0.0154419| var\_p    |   0.5825667| Peak\_m          |  -0.4191859|
| jsd\_p    |  -0.1826259| mae\_m    |   0.5186601| Peak\_p          |  -0.2512678|
| kl\_m     |   0.5504755| mae\_p    |   0.5321362| ImpulseFactor\_m |  -0.5075726|
| kl\_p     |   0.6731723| rmse\_m   |   0.5857615| ImpulseFactor\_p |  -0.4730530|
| arclen\_m |  -0.0765462| rmse\_p   |   0.6123484| NA               |          NA|

The correlation matrix looks as in figure .

<img src="fire-drill-technical-report_files/figure-markdown_github/p2-no-corr-mat-1.png" alt="Correlation matrix for non-aligned scores using pattern II."  />
<p class="caption">
Correlation matrix for non-aligned scores using pattern II.
</p>

While some correlations are weaker now, we observe more agreement
between the scores, i.e., more correlations tend to be positive. *Peak*
and *Impulse-factor* however have negative correlations now.

Pattern III (average)
---------------------

The 3rd pattern that was produced as weighted average over all ground
truth is scored in this section. **Note!** There is one important
difference here: the weighted-average pattern does not have the same
intervals as our initial pattern – in fact, we cannot make any
assumptions about any of the intervals. Therefore, we will compute this
align with **ten equally-long** intervals. This amount was chosen
arbitrarily as a trade-off between the time it takes to compute, and the
resolution of the results. Adding more intervals increases both,
computing time and resolution exponentially, however the latter much
less.

``` r
library(foreach)

p3_avg_align <- loadResultsOrCompute(file = "../results/p3_avg_align.rds", computeExpr = {
  # Let's compute all projects in parallel!
  cl <- parallel::makePSOCKcluster(length(project_signals))
  unlist(doWithParallelClusterExplicit(cl = cl, expr = {
    foreach::foreach(
      projectName = names(project_signals),
      .inorder = FALSE,
      .packages = c("parallel")
    ) %dopar% {
      source("./common-funcs.R")
      source("../models/modelsR6.R")
      source("../models/SRBTW-R6.R")
      
      cl_nested <- parallel::makePSOCKcluster(5)
      `names<-`(list(doWithParallelClusterExplicit(cl = cl_nested, expr = {
        temp <- time_warp_project(
          thetaB = seq(from = 0, to = 1, by = 0.1), # important!
          pattern = p3_avg_signals, project = project_signals[[projectName]])
        temp$fit(verbose = TRUE)
        temp # return the instance, it includes the FitResult
      })), projectName)
    }
  }))
})
```

``` r
p3_avg_scores <- loadResultsOrCompute(file = "../results/p3_avg_scores.rds", computeExpr = {
  as.data.frame(compute_all_scores(alignment = p3_avg_align, patternName = "p3_avg"))
})
```

|                  |     pr\_1|     pr\_2|   pr\_3|     pr\_4|        pr\_5|          pr\_6|        pr\_7|    pr\_8|        pr\_9|
|:-----------------|---------:|---------:|-------:|---------:|------------:|--------------:|------------:|--------:|------------:|
| area\_m          |      0.93|      0.90|    0.88|      0.94|         0.94|           0.93|         0.89|     0.93|         0.92|
| area\_p          |      0.73|      0.66|    0.60|      0.77|         0.79|           0.74|         0.62|     0.75|         0.71|
| corr\_m          |      0.61|      0.65|    0.53|      0.73|         0.76|           0.86|         0.70|     0.78|         0.73|
| corr\_p          |      0.13|      0.14|    0.07|      0.28|         0.32|           0.55|         0.22|     0.37|         0.27|
| jsd\_m           |      0.45|      0.43|    0.38|      0.56|         0.52|           0.54|         0.39|     0.53|         0.48|
| jsd\_p           |      0.04|      0.03|    0.02|      0.08|         0.07|           0.06|         0.02|     0.07|         0.04|
| kl\_m            |     46.44|     23.01|    4.43|     38.29|      1023.40|        4098.55|      6731.77|    10.23|      2921.25|
| kl\_p            |  28603.90|  10875.90|  190.61|  29987.58|  12752649.76|  323943591\.52|  57623738.67|  1627.21|  26054558.68|
| arclen\_m        |      0.57|      0.68|    0.53|      0.66|         0.63|           0.65|         0.80|     0.74|         0.60|
| arclen\_p        |      0.10|      0.19|    0.07|      0.17|         0.12|           0.14|         0.31|     0.26|         0.13|
| sd\_m            |      0.83|      0.84|    0.79|      0.86|         0.89|           0.88|         0.83|     0.87|         0.86|
| sd\_p            |      0.48|      0.50|    0.38|      0.54|         0.61|           0.60|         0.47|     0.58|         0.55|
| var\_m           |      0.97|      0.97|    0.95|      0.98|         0.98|           0.98|         0.97|     0.98|         0.98|
| var\_p           |      0.88|      0.90|    0.81|      0.90|         0.94|           0.93|         0.88|     0.93|         0.91|
| mae\_m           |      0.93|      0.90|    0.88|      0.94|         0.94|           0.93|         0.89|     0.93|         0.92|
| mae\_p           |      0.73|      0.66|    0.60|      0.77|         0.79|           0.74|         0.62|     0.75|         0.71|
| rmse\_m          |      0.88|      0.87|    0.84|      0.90|         0.92|           0.91|         0.85|     0.90|         0.89|
| rmse\_p          |      0.59|      0.58|    0.50|      0.65|         0.70|           0.67|         0.51|     0.67|         0.62|
| RMS\_m           |      0.54|      0.50|    0.48|      0.67|         0.64|           0.55|         0.50|     0.58|         0.62|
| RMS\_p           |      0.07|      0.05|    0.04|      0.19|         0.16|           0.09|         0.04|     0.09|         0.13|
| Kurtosis\_m      |      0.13|      0.37|    0.18|      0.28|         0.26|           0.16|         0.20|     0.43|         0.19|
| Kurtosis\_p      |      0.00|      0.00|    0.00|      0.00|         0.00|           0.00|         0.00|     0.01|         0.00|
| Peak\_m          |      0.53|      0.71|    0.57|      0.60|         0.68|           0.58|         0.47|     0.72|         0.54|
| Peak\_p          |      0.07|      0.20|    0.10|      0.13|         0.20|           0.11|         0.03|     0.21|         0.08|
| ImpulseFactor\_m |      0.61|      0.74|    0.69|      0.62|         0.67|           0.44|         0.49|     0.70|         0.42|
| ImpulseFactor\_p |      0.12|      0.28|    0.21|      0.14|         0.18|           0.03|         0.02|     0.22|         0.03|

The correlation of just the ground truth with all scores is in table .

| Score     |       Value| Score     |       Value| Score            |       Value|
|:----------|-----------:|:----------|-----------:|:-----------------|-----------:|
| area\_m   |  -0.1308766| arclen\_p |  -0.2487774| RMS\_m           |   0.2993487|
| area\_p   |  -0.1294226| sd\_m     |  -0.3115364| RMS\_p           |   0.3621570|
| corr\_m   |  -0.2236643| sd\_p     |  -0.3275009| Kurtosis\_m      |  -0.3521132|
| corr\_p   |  -0.1652319| var\_m    |  -0.4254858| Kurtosis\_p      |  -0.4821234|
| jsd\_m    |   0.0300469| var\_p    |  -0.4289845| Peak\_m          |  -0.4520734|
| jsd\_p    |  -0.0243270| mae\_m    |  -0.1309105| Peak\_p          |  -0.4502977|
| kl\_m     |   0.0246984| mae\_p    |  -0.1296077| ImpulseFactor\_m |  -0.2753840|
| kl\_p     |  -0.1039983| rmse\_m   |  -0.2254399| ImpulseFactor\_p |  -0.2790409|
| arclen\_m |  -0.2946963| rmse\_p   |  -0.2372254| NA               |          NA|

The correlation matrix looks as in figure .

<img src="fire-drill-technical-report_files/figure-markdown_github/p3-avg-corr-mat-1.png" alt="Correlation matrix for scores using pattern III (average)."  />
<p class="caption">
Correlation matrix for scores using pattern III (average).
</p>

I suppose that the most significant result here is the positive
Jensen–Shannon divergence score correlation.

Pattern III (average, no alignment)
-----------------------------------

Before we go any further, I would also like to compute the scores for
this data-driven pattern **without** having the projects aligned. After
manually inspecting some of these alignments, it turns out that some are
quite extreme. Since this pattern is a weighted average over all ground
truth, not much alignment should be required.

``` r
# We'll have to mimic an 'aligned' object, which is a list
# of srBTAW_MultiVartype instances. We can clone it and just
# undo the time warping.
p3_avg_no_align <- list()

for (project in ground_truth$project) {
  temp <- p3_avg_align[[project]]$clone()
  temp$setParams(params = `names<-`(rep(1 / temp$getNumParams(),
                  temp$getNumParams()), temp$getParamNames()))
  p3_avg_no_align[[project]] <- temp
}
```

``` r
p3_avg_no_scores <- loadResultsOrCompute(file = "../results/p3_avg_no_scores.rds", computeExpr = {
  as.data.frame(compute_all_scores(alignment = p3_avg_no_align, patternName = "p3_avg"))
})
```

|                  |   pr\_1|    pr\_2|     pr\_3|         pr\_4|    pr\_5|        pr\_6|      pr\_7|   pr\_8|     pr\_9|
|:-----------------|-------:|--------:|---------:|-------------:|--------:|------------:|----------:|-------:|---------:|
| area\_m          |    0.80|     0.88|      0.91|  9.000000e-01|     0.86|         0.88|       0.91|    0.84|      0.91|
| area\_p          |    0.39|     0.59|      0.69|  6.600000e-01|     0.55|         0.60|       0.67|    0.50|      0.69|
| corr\_m          |    0.66|     0.50|      0.71|  8.200000e-01|     0.49|         0.58|       0.73|    0.59|      0.77|
| corr\_p          |    0.19|     0.06|      0.22|  4.400000e-01|     0.05|         0.11|       0.27|    0.08|      0.34|
| jsd\_m           |    0.36|     0.37|      0.49|  4.900000e-01|     0.33|         0.36|       0.46|    0.39|      0.46|
| jsd\_p           |    0.01|     0.01|      0.04|  4.000000e-02|     0.01|         0.01|       0.03|    0.02|      0.03|
| kl\_m            |    7.33|    17.22|     60.96|  5.231730e+03|     8.73|       290.40|      95.40|    4.92|     19.84|
| kl\_p            |  604.15|  7730.05|  57459.29|  6.591741e+09|  1016.03|  11895835.33|  105525.86|  386.95|  18235.65|
| arclen\_m        |    0.75|     0.61|      0.71|  6.700000e-01|     0.58|         0.63|       0.65|    0.71|      0.76|
| arclen\_p        |    0.26|     0.13|      0.23|  1.600000e-01|     0.10|         0.15|       0.15|    0.21|      0.29|
| sd\_m            |    0.75|     0.77|      0.83|  8.700000e-01|     0.75|         0.80|       0.84|    0.80|      0.84|
| sd\_p            |    0.30|     0.34|      0.45|  5.600000e-01|     0.31|         0.40|       0.48|    0.40|      0.50|
| var\_m           |    0.93|     0.94|      0.96|  9.800000e-01|     0.93|         0.95|       0.96|    0.95|      0.97|
| var\_p           |    0.73|     0.77|      0.85|  9.100000e-01|     0.74|         0.82|       0.86|    0.82|      0.89|
| mae\_m           |    0.80|     0.88|      0.91|  9.000000e-01|     0.86|         0.88|       0.91|    0.84|      0.91|
| mae\_p           |    0.39|     0.59|      0.69|  6.600000e-01|     0.55|         0.60|       0.67|    0.50|      0.69|
| rmse\_m          |    0.73|     0.82|      0.87|  8.900000e-01|     0.80|         0.84|       0.88|    0.80|      0.88|
| rmse\_p          |    0.27|     0.45|      0.57|  6.100000e-01|     0.41|         0.50|       0.58|    0.40|      0.60|
| RMS\_m           |    0.46|     0.43|      0.50|  5.500000e-01|     0.39|         0.48|       0.56|    0.47|      0.54|
| RMS\_p           |    0.02|     0.03|      0.06|  9.000000e-02|     0.02|         0.05|       0.09|    0.03|      0.07|
| Kurtosis\_m      |    0.14|     0.11|      0.21|  1.200000e-01|     0.12|         0.10|       0.32|    0.37|      0.37|
| Kurtosis\_p      |    0.00|     0.00|      0.00|  0.000000e+00|     0.00|         0.00|       0.00|    0.00|      0.00|
| Peak\_m          |    0.47|     0.53|      0.60|  5.800000e-01|     0.57|         0.54|       0.68|    0.71|      0.72|
| Peak\_p          |    0.03|     0.07|      0.13|  1.100000e-01|     0.10|         0.08|       0.20|    0.20|      0.21|
| ImpulseFactor\_m |    0.49|     0.71|      0.69|  3.700000e-01|     0.78|         0.37|       0.74|    0.76|      0.69|
| ImpulseFactor\_p |    0.05|     0.17|      0.22|  1.000000e-02|     0.34|         0.01|       0.29|    0.30|      0.20|

The correlation of just the ground truth with all scores is in table .

| Score     |      Value| Score     |      Value| Score            |       Value|
|:----------|----------:|:----------|----------:|:-----------------|-----------:|
| area\_m   |  0.6670842| arclen\_p |  0.2514802| RMS\_m           |   0.7467092|
| area\_p   |  0.6894208| sd\_m     |  0.8218241| RMS\_p           |   0.7740284|
| corr\_m   |  0.8412487| sd\_p     |  0.8436768| Kurtosis\_m      |   0.0543283|
| corr\_p   |  0.8835587| var\_m    |  0.7916103| Kurtosis\_p      |   0.0289034|
| jsd\_m    |  0.8823326| var\_p    |  0.8025434| Peak\_m          |   0.1879883|
| jsd\_p    |  0.8724026| mae\_m    |  0.6670619| Peak\_p          |   0.1986320|
| kl\_m     |  0.6774499| mae\_p    |  0.6893942| ImpulseFactor\_m |  -0.3763130|
| kl\_p     |  0.6729326| rmse\_m   |  0.7436967| ImpulseFactor\_p |  -0.3217380|
| arclen\_m |  0.2914946| rmse\_p   |  0.7784040| NA               |          NA|

The correlation matrix looks as in figure .

<img src="fire-drill-technical-report_files/figure-markdown_github/p3-avg-no-corr-mat-1.png" alt="Correlation matrix for non-aligned scores using pattern III (average)."  />
<p class="caption">
Correlation matrix for non-aligned scores using pattern III (average).
</p>

And here in figure we got the result that was the most-expected. We get
almost always positive correlations, and most of them range between
medium and significant strength. If we look at the correlation for
`sd_p`, it is almost perfect with 0.844. The Jensen–Shannon divergence
score (note: while called divergence, this is a score, the higher the
better) of 0.882 is at a high level now. I mention this because this
measure is less primitive than the others and tends to capture more
properties of the signals. If this is high, it means we get a robust
measure that ought to be usable stand-alone. The other low-level scores
probably would need to be combined instead.

### Linear combination of scores

So far we have tested whether the calculated scores correlate with the
scores of the ground truth, and we find some good examples. However,
these scores are not scaled or translated in any way, so it is probably
best to A) take multiple scores into account and B) create a regression
model that makes these adjustments. We will test some linear combination
of the scores `corr_p`, `jsd_m`, `RMS_m` and `sd_p`.

``` r
p3_avg_lm <- stats::glm(formula = gt_consensus ~ corr_p + jsd_m + RMS_m + sd_p, data = temp)
stats::coef(p3_avg_lm)
```

    ## (Intercept)      corr_p       jsd_m       RMS_m        sd_p 
    ##   0.1127856   1.4429241   2.3439285  -3.2697287   1.2561099

``` r
plot(p3_avg_lm, ask = FALSE, which = 1:2)
```

![](fire-drill-technical-report_files/figure-markdown_github/unnamed-chunk-56-1.png)![](fire-drill-technical-report_files/figure-markdown_github/unnamed-chunk-56-2.png)

Of course, this model should not be used to predict on the training
data, but what we wanted to learn here is simply how to linearly combine
the scores in order to get scores that are in the same range as the
ground truth, we learn a re-scale so to say.

``` r
p3_avg_lm_scores <- stats::predict(p3_avg_lm, temp)
round(p3_avg_lm_scores * 10, 3)
```

    ## project_1 project_2 project_3 project_4 project_5 project_6 project_7 project_8 
    ##     1.087     0.827     4.968     7.975     0.739     0.311     3.565     1.180 
    ## project_9 
    ##     5.348

``` r
stats::cor(p3_avg_lm_scores, ground_truth$consensus_score)
```

    ## [1] 0.9485156

This also increased the correlation to 0.949.

Pattern III (b)
---------------

The third and last pattern is based on the ground truth only. Starting
with straight lines and equally long intervals, a pattern was generated
and selected using a best trade-off between number of parameters and
highest *likelihood*, or by the lowest loss. Like pattern II, all of
these patterns were produced with time warping applied, so that we do
not need to align the projects to it. We have produced 16 patterns, so
let’s compute the scores for all of them.

``` r
p3b_no_scores <- loadResultsOrCompute(file = "../results/p3b_no_scores.rds", computeExpr = {
  unlist(doWithParallelCluster(numCores = min(8, parallel::detectCores()), expr = {
    foreach::foreach(
      numIntervals = seq_len(length.out = 16),
      .inorder = FALSE
    ) %dopar% {
      source("./common-funcs.R")
      source("../models/modelsR6.R")
      source("../models/SRBTW-R6.R")
      
      temp_no_align <- list()
      
      for (project in ground_truth$project) {
        inst <- srBTAW$new(
          theta_b = seq(0, 1, length.out = 2),
          gamma_bed = c(0, 1, .Machine$double.eps),
          lambda = .Machine$double.eps,
          begin = 0, end = 1, openBegin = FALSE, openEnd = FALSE,
          useAmplitudeWarping = FALSE)
        
        for (vartype in names(weight_vartype)) {
          # Set WP and WC:
          wp <- p3b_all[[paste0("i_", numIntervals)]]$signals[[vartype]]
          wc <- project_signals[[project]][[vartype]]
            
          inst$setSignal(signal = wp)
          inst$setSignal(signal = wc)
          
          srbtwbaw <- inst$.__enclos_env__$private$createInstance(
            wp = wp$get0Function(), wc = wc$get0Function())
          inst$.__enclos_env__$private$addInstance(
            instance = srbtwbaw, wpName = wp$getName(), wcName = wc$getName())
        }
        
        inst$setParams(params = c(vtl_1 = 1))
        temp_no_align[[project]] <- inst
      }
      
      res <- `names<-`(list(tryCatch(expr = {
        # compute_all_scores is already parallel!
        temp <-as.data.frame(compute_all_scores(alignment = temp_no_align, patternName = "p3"))
        saveRDS(object = temp, file = paste0("../results/p3-compute/p3b_no_", numIntervals, ".rds"))
        temp
      }, error = function(cond) {
        list(cond = cond, tb = traceback())
      })), paste0("i_", numIntervals))
    }
  }), recursive = FALSE)
})
```

We will have a lot of results. In order to give an overview, we will
show the correlation of the ground truth with the scores as obtained by
scoring against each of the 16 patterns. The correlation plot in figure
should give us a good idea of how the scores changes with increasing
number of parameters.

<img src="fire-drill-technical-report_files/figure-markdown_github/p3b-no-scores-corr-1.png" alt="Correlation-scores for all 16 patterns (type III, b)."  />
<p class="caption">
Correlation-scores for all 16 patterns (type III, b).
</p>

If `RMS` was to score to use, then the even the 1-interval pattern will
do, as the correlation for it is always strongly positive. In general,
we can observe how some correlations get weaker, and some get stronger
for patterns with higher number of parameters. There is no clear winner
here, and the results are quite similar. What can say clearly is, that
it is likely not worth to use highly parameterized models to detect the
Fire Drill as it was present in our projects, as the manifestation is
just not strong enough to warrant for patterns with high degrees of
freedom. It is probably best, to use one of the other pattern types.

Let’s also show an overview of the correlation with the ground truth for
all of the other patterns:

<img src="fire-drill-technical-report_files/figure-markdown_github/pall-corr-1.png" alt="Overview of correlation-scores for all other types of patterns."  />
<p class="caption">
Overview of correlation-scores for all other types of patterns.
</p>

The correlation overview in figure suggests that the no-alignment
patterns have most of the strong positive-correlated scores. However, it
appears that it was still worth adapting our initial pattern using the
ground truth, as it is the only pattern with very high positive
correlations for Peak and Impulse-factor. These however disappear, if we
do not do the double warping.

### Linear combination of scores

The pattern that incurred the lowest loss was number **4**. We should
however also test the pattern that had the highest correlations
(positive *or* negative) on average:

``` r
p3b_corr_all <- apply(X = p3b_corr, MARGIN = 1,
                      FUN = function(row) mean(abs(row)))
```

|                         |       corr|
|:------------------------|----------:|
| P III (b) \[numInt=1\]  |  0.3303488|
| P III (b) \[numInt=2\]  |  0.3960491|
| P III (b) \[numInt=3\]  |  0.3836158|
| P III (b) \[numInt=4\]  |  0.3516578|
| P III (b) \[numInt=5\]  |  0.3421378|
| P III (b) \[numInt=6\]  |  0.3365299|
| P III (b) \[numInt=7\]  |  0.3421056|
| P III (b) \[numInt=9\]  |  0.3440612|
| P III (b) \[numInt=10\] |  0.2858276|
| P III (b) \[numInt=11\] |  0.4983041|
| P III (b) \[numInt=12\] |  0.2555541|
| P III (b) \[numInt=13\] |  0.4151310|
| P III (b) \[numInt=14\] |  0.3560980|
| P III (b) \[numInt=15\] |  0.4242694|
| P III (b) \[numInt=16\] |  0.3373230|

In table we show the mean absolute correlation for the scores of all
projects as computed against each pattern of type III (b).

``` r
p3b_highest_corr <- names(which.max(p3b_corr_all))
p3b_highest_corr
```

    ## [1] "P III (b) [numInt=11]"

The pattern using **11** has the highest correlation.

Taking the best pattern of type III (b), we can also do a linear
combination again. “Best” may refer to the pattern with the highest
correlations over all scores or the one with the lowest loss when
fitting to all projects.

``` r
# Take the project with the highest mean absolute correlation.
temp <- cbind(
  data.frame(gt_consensus = ground_truth$consensus_score),
  p3b_no_scores[[paste0("i_", gsub("[^0-9]", "", p3b_highest_corr))]])
p3b_lm <- stats::lm(
  formula = gt_consensus ~ area_p + corr_p + jsd_m + sd_p + rmse_p + RMS_p,
  data = temp)
stats::coef(p3b_lm)
```

    ## (Intercept)      area_p      corr_p       jsd_m        sd_p      rmse_p 
    ##    2.616642  -38.644714    3.048282   12.414117  -18.992103   50.599212 
    ##       RMS_p 
    ##  -46.685270

``` r
plot(p3b_lm, ask = FALSE, which = 1:2)
```

![](fire-drill-technical-report_files/figure-markdown_github/unnamed-chunk-61-1.png)![](fire-drill-technical-report_files/figure-markdown_github/unnamed-chunk-61-2.png)

``` r
p3b_lm_scores <- stats::predict(p3b_lm, temp)
round(p3b_lm_scores * 10, 3)
```

    ## project_1 project_2 project_3 project_4 project_5 project_6 project_7 project_8 
    ##     0.941    -0.573     6.477     7.214     0.978     2.972     3.473     0.006 
    ## project_9 
    ##     4.511

``` r
stats::cor(p3b_lm_scores, ground_truth$consensus_score)
```

    ## [1] 0.9798728

With this linear model, we can report a high correlation with the
consensus, too.

References
==========

Akaike, Hirotugu. 1981. “Likelihood of a Model and Information
Criteria.” *Journal of Econometrics* 16 (1): 3–14.

Brown, William J, Raphael C Malveau, Hays W McCormick III, and Thomas J
Mowbray. 1998. “Refactoring Software, Architectures, and Projects in
Crisis.” John Wiley; Sons, Inc, Canada.

Cohen, Jacob. 1968. “Weighted Kappa: Nominal Scale Agreement Provision
for Scaled Disagreement or Partial Credit.” *Psychological Bulletin* 70
(4): 213.

Hönel, Sebastian. 2020. “Git Density 2020.2: Analyze Git Repositories to
Extract the Source Code Density and Other Commit Properties,” February.
<https://doi.org/10.5281/zenodo.3662768>.

Hönel, Sebastian, Morgan Ericsson, Welf Löwe, and Anna Wingkvist. 2020.
“Using Source Code Density to Improve the Accuracy of Automatic Commit
Classification into Maintenance Activities.” *Journal of Systems and
Software*, 110673.

Landis, J Richard, and Gary G Koch. 1977. “An Application of
Hierarchical Kappa-Type Statistics in the Assessment of Majority
Agreement Among Multiple Observers.” *Biometrics*, 363–74.

Xi, Fengfeng, Qiao Sun, and Govindappa Krishnappa. 2000. “Bearing
Diagnostics Based on Pattern Recognition of Statistical Parameters.”
*Journal of Vibration and Control* 6 (3): 375–92.
<http://www.acoustics.asn.au/conference_proceedings/ICSVS-1997/pdf/scan/sv970356.pdf>.

[1] <a href="https://github.com/sse-lnu/anti-pattern-models/blob/master/notebooks/comm-class-models.Rmd" class="uri">https://github.com/sse-lnu/anti-pattern-models/blob/master/notebooks/comm-class-models.Rmd</a>
