# Detecting the Fire Drill anti-pattern using Source Code and issue-tracking data

In this repository, we develop methods that can model and detect the
presence of so-called anti-patterns (AP). In this repository, you will
find:

-   [**Models** (R)](./models) - Multivariate models (arbitrary many)
    that can encapsulate the definition of an anti-pattern, using
    arbitrary user-defined intervals and losses. The models can then be
    fit to the data (or vice versa). Also, there is lots of
    functionality for **quantifying differences**, esp. for **scoring**.
-   [**Notebooks** (R)](./notebooks) - Notebooks that can be (re-)run by
    users to reproduce our results. All results are included, such that
    the notebooks will only recompute them if you delete them. The
    notebooks are very detailed and document all steps necessary. See
    (Hönel 2022) for a pre-rendered PDF.
-   [**Data** (CSV)](./data) and (precomputed) [**Results**
    (RDS)](./results) - All data required for reproduction is included.
    All the results, too. Some of them take days to compute, so be
    aware. Also see (Hönel et al. 2022).

There is a case study that makes extensive use of the data and models
(Picha et al. 2022).

# Abstract

Detecting the presence of project management anti-patterns (AP)
currently requires experts on the matter and is an expensive endeavor.
Worse, experts may introduce their individual subjectivity or bias.
Using the Fire Drill AP, we first introduce a novel way to translate
descriptions into detectable AP that are comprised of arbitrary metrics
and events such as logged time or maintenance activities, which are
mined from the underlying source code or issue-tracking data, thus
making the description objective as it becomes data-based. Secondly, we
demonstrate a novel method to quantify and score the deviations of
real-world projects to data-based AP descriptions. Using fifteen
real-world projects that exhibit a Fire Drill to some degree, we show
how to further enhance the translated AP. The ground truth in these
projects was extracted from two individual experts and consensus was
found between them. Our evaluation spans four kinds of patterns, where
the first is purely derived from description, the second type is
enhanced by data, and the third kind is derived from data only. The
fourth type then is a derivative meta-process pattern. We introduce a
novel method called automatic calibration, that optimizes a pattern such
that only necessary and important scores remain that suffice to
confidently detect the degree to which the AP is present. Without
automatic calibration, the proposed patterns show only weak potential
for detecting the presence. Enriching the AP with data from real-world
projects significantly improves the potential. We conclude that the
presence of similar patterns is most certainly detectable. Furthermore,
any pattern that can be characteristically modeled using the proposed
approach is potentially well detectable.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-honel2021technical" class="csl-entry">

Hönel, Sebastian. 2022. “Technical Reports Compilation: Detecting the
Fire Drill Anti-Pattern Using Source Code and Issue-Tracking Data.”
*CoRR*. <https://doi.org/10.48550/arXiv.2104.15090>.

</div>

<div id="ref-honel_picha_2021" class="csl-entry">

Hönel, Sebastian, Petr Pícha, Premek Brada, and Lenka Rychtarova. 2022.
“Detection of the Fire Drill Anti-Pattern: 15 Real-World Projects with
Ground Truth, Issue-Tracking Data, Source Code Density, Models and
Code.” Zenodo. <https://doi.org/10.5281/zenodo.5992621>.

</div>

<div id="ref-picha2022Firedrill" class="csl-entry">

Picha, Petr, Sebastian Hönel, Premek Brada, Morgan Ericsson, Welf Löwe,
Anna wingkvist, and Jakub Danek. 2022. “Process Anti-Pattern Detection
in Student Projects – a Case Study.” In *EuroPLoP’22: European
Conference on Pattern Languages of Programs 2022, Irsee, Deutschland,
July 6 - 10, 2022*. ACM.

</div>

</div>
