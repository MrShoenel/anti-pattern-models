# Detecting the Fire Drill Anti-pattern Using Source Code and Issue-Tracking Data [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6884866.svg)](https://doi.org/10.5281/zenodo.6884866)

In this repository, we develop methods that can model and detect the
presence of so-called anti-patterns (AP). In this repository, you will
find:

- [**Models** (R)](./models) - Multivariate models (arbitrary many) that
  can encapsulate the definition of an anti-pattern, using arbitrary
  user-defined intervals and losses. The models can then be fit to the
  data (or vice versa). Also, there is lots of functionality for
  **quantifying differences**, esp. for **scoring**.
- [**Notebooks** (R)](./notebooks) - Notebooks that can be (re-)run by
  users to reproduce our results. All results are included, such that
  the notebooks will only recompute them if you delete them. The
  notebooks are very detailed and document all steps necessary. See
  (Hönel 2023b) for a pre-rendered PDF.
- [**Data** (CSV)](./data) and (precomputed) [**Results**
  (RDS)](./results) - All data required for reproduction is included.
  All the results, too. Some of them take days to compute, so be aware.
  Also see (Hönel, Pícha, et al. 2023).

There is a pilot study that makes extensive use of the data and models
(Picha et al. 2022). The eight version of the repository, dataset, and
technical report serve then as the basis for a proper, subsequent
embedded case study (Hönel, Picha, et al. 2023). This repository has a
release on Zenodo (Hönel 2023a). From version seven and onwards, each
technical reports compilation will be paired with a separate release of
this repository on Zenodo.

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
found between them. We introduce a novel method called automatic
calibration, that optimizes a pattern such that only necessary and
important scores remain that suffice to confidently detect the degree to
which the AP is present. Without automatic calibration, the proposed
patterns show only weak potential for detecting the presence. Enriching
the AP with data from real-world projects significantly improves the
potential. We also introduce a no-pattern approach that exploits the
ground truth for establishing a new, quantitative understanding of the
phenomenon, as well as for finding gray-/black-box predictive models. We
conclude that the presence detection and severity assessment of the Fire
Drill anti-pattern, as well as some of its related and similar patterns,
is certainly possible using some of the presented approaches.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-gitHub_repo_latest" class="csl-entry">

Hönel, Sebastian. 2023a. “<span class="nocase">GitHub Repository:
Detecting the Fire Drill Anti-pattern Using Source Code and
Issue-Tracking Data</span>,” January.
<https://doi.org/10.5281/zenodo.6884866>.

</div>

<div id="ref-honel2021technical" class="csl-entry">

———. 2023b. “Technical Reports Compilation: Detecting the Fire Drill
Anti-Pattern Using Source Code and Issue-Tracking Data.” *CoRR*. arXiv.
<https://doi.org/10.48550/arXiv.2104.15090>.

</div>

<div id="ref-honel2023embedded" class="csl-entry">

Hönel, Sebastian, Petr Picha, Morgan Ericsson, Premek Brada, Welf Löwe,
and Anna Wingkvist. 2023. “Activity-Based Detection of Pattern-Like
Phenomena: An Embedded Case Study of the Fire Drill.” *E-Informatica
Software Engineering Journal*, February.

</div>

<div id="ref-honel_picha_2021" class="csl-entry">

Hönel, Sebastian, Petr Pícha, Premek Brada, Lenka Rychtarova, and Jakub
Danek. 2023. “<span class="nocase">Detection of the Fire Drill
anti-pattern: 15 real-world projects with ground truth, issue-tracking
data, source code density, models and code</span>.” Zenodo.
<https://doi.org/10.5281/zenodo.4734053>.

</div>

<div id="ref-picha2022Firedrill" class="csl-entry">

Picha, Petr, Sebastian Hönel, Premek Brada, Morgan Ericsson, Welf Löwe,
Anna Wingkvist, and Jakub Danek. 2022. “Process Anti-Pattern Detection
in Student Projects – a Case Study.” In *EuroPLoP’22: European
Conference on Pattern Languages of Programs 2022, Irsee, Deutschland,
July 6 - 10, 2022*. ACM. <https://doi.org/10.1145/3551902.3551965>.

</div>

</div>
