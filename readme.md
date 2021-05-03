Detecting the Fire Drill anti-pattern using Source Code
=======================================================

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
    (Hönel 2021) for a pre-rendered PDF.
-   [**Data** (CSV)](./data) and (precomputed) [**Results**
    (RDS)](./results) - All data required for reproduction is included.
    All the results, too. Some of them take days to compute, so be
    aware. Also see (Hönel et al. 2021).

Abstract
========

Detecting the presence of project management anti-patterns (AP)
currently requires experts on the matter and is an expensive endeavor.
Worse, experts may introduce their individual subjectivity or bias.
Using the Fire Drill AP, we first introduce a novel way to translate
descriptions into detectable AP that are comprised of arbitrary metrics
and events such as maintenance activities, which are mined from the
underlying source code, thus making the description objective as it
becomes data-based. Secondly, we demonstrate a novel method to quantify
and score the deviations of real-world projects to data-based AP
descriptions. Using nine real-world projects that exhibit a Fire Drill
to some degree, we show how to further enhance the translated AP. The
ground truth in these projects was extracted from two individual experts
and consensus was found between them. Our evaluation spans three kinds
of pattern, where the first is purely derived from description, the
second type is enhanced by data, and the third kind is derived from data
only. The Fire Drill AP as translated from description only shows weak
potential of confidently detecting the presence of the anti-pattern in a
project. Enriching the AP with data from real-world projects
significantly improves the detection. Using patterns derived from data
only leads to almost perfect correlations of the scores with the ground
truth. Some APs share symptoms with the Fire Drill AP, and we conclude
that the presence of similar patterns is most certainly detectable.
Furthermore, any pattern that can be characteristically modelled using
the proposed approach is potentially well detectable.

References
==========

Hönel, Sebastian. 2021. “Technical Reports Compilation: Detecting the
Fire Drill Anti-Pattern Using Source Code.”
<http://arxiv.org/abs/2104.15090>.

Hönel, Sebastian, Petr Pícha, Premek Brada, and Lenka Rychtarova. 2021.
“Detection of the Fire Drill Anti-Pattern: Nine Real-World Projects with
Ground Truth, Issue-Tracking Data, Source Code Density, Models and
Code.” Zenodo. <https://doi.org/10.5281/zenodo.4734053>.
