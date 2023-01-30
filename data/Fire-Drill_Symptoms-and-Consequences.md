---
title: "Fire Drill: The Anti-pattern"
author: "Sebastian Hönel, Petr Picha, Premek Brada"
geometry: margin=2.5cm
papersize: A4
fontsize: 11pt
documentclass: scrartcl
classoption: "egregdoesnotlikesansseriftitles"
numbersections: true
subparagraph: yes # see https://bit.ly/2rzAIgv
colorlinks: true # Set to true and change the 3 colors below if you like
linkcolor: "Maroon" # Default is "Maroon"
citecolor: "Blue" # Default is "Blue"
urlcolor: "Blue" # Default is "Blue"
linestretch: 1.1 # For better readability
biblio-title: "References" # Set title of the References/Bibliography
biblatexoptions: "backend=biber,style=numeric,natbib=true"
---


**Note**: This document has a canonical URL: <https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/catalogue/Fire_Drill.md>.

________

[Home](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/README.md) > [Catalogue](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/Antipatterns_catalogue.md) > Fire Drill


# Fire Drill

## Also Known As

n/a

## Summary (partly solution)

Requirements and Analysis phases prolonged and consuming disproportionate amount of resources (because management want to do them "right"), then frantic "everything needs to be done yesterday" period to finish on time (when management finds out they wasted most of project's schedule and resources on analysis).


## Solution

With only a fraction of the required resources (e.g., time) left,
- `S1`. Pressure the team to meet the original deadline(s), and
- `S2`. Pressure the team to ship the full, initially agreed-upon deliverables.


## Context

Waterfall(ish) projects, especially when project oversight is loose and/or management is not driven by outcome.

## Unbalanced Forces

- `F1`. need (desire) to have specifications perfect
- `F2`. management consumed by internal (political) issues
- `F3`. actual development of a high-quality product takes time
- `F4`. quality objectives formally stated and high
- `F5`. strict deadlines for delivery

## Symptoms and Consequences

For each symptom/consequence, we also add a sublist with empirically gathered, concrete instances (marked by `E[0-9]+`). Only `SC1` through `SC7` are the original symptoms and consequences, everything else was added manually.
The empirically found concrete instances were found when analyzing the data [[HON'23]](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/References.md).

The empirical problem instances are ascribed a __severity__ using the following Likert scale:

- `[0]` __None__: not at all a problem: applies mostly to false positives (e.g., a typical symptom that was caused out of the studied context and had no adverse effects)
- `[1]` __Miniscule__: only slight indications of typical symptom(s) identified by at least one rater
- `[2]` __Minor__: multiple indicators and/or measurable/documented symptom(s); sometimes corroborated by another rater
- `[3]` __Moderate__: clearly identifiable and reoccurring symptoms or direct corroboration
- `[4]` __Significant__: like moderate, but the higher severity is evident through additional data- or observer-triangulation
- `[5]` __Serious__: agreement on the (recurrent) severe presence of a symptom by observer-triangulation (often all raters) and/or data-triangulation

The severity attached follows the format `[project, rater(s), severity]`, e.g., `[1,ab,3]` indicates that raters A and B identify a problem instance in project `1`, and the severity is `3/5`.

The empirical observations were only logged if the rater's notes allowed it. For example, some notes indicate a problem, without a cause: _"descope in later stages of the project"_ or _"poor testing"_. Here, we cannot assign an instance of symptom/cause.
Also, if the severity cannot be decided between two levels, then the rater's Fire Drill severity is applied to make a decision.

Finally, some observations are aggregated by observer triangulation (using the expression `[..], [..], .. -> [..]`).
The general rule is to upgrade the severity by `+1` if at least two observers find the same problem instance in a project. If all three observers agree, the severity is set to `5`. Here are some examples: `1+1 -> 2, 1+2 -> 3, 1+3 -> 3, 2+3 -> 4, 3+4 -> 5, 3+3+4 -> 5`.

- `SC1`. long period at project start where activities connected to requirements, analysis and planning prevail, and design and implementation activities are rare
  - `E01`. can be caused by, e.g., the team needing some time for familiarization with the (new, changed) tools, way of communication, or process in the beginning `[1,A,1]`, `[3,B,2]`, `[5,A,1]`, `[7,B,1]`, `[9,B,2]`, `[13,A,1]`, `[14,A,1]`, `[15,A,2]`
  - `E02`. work delayed due to external factors, such as ramping-down of previous project or other, non-related work `[1,A,0]`, `[15,A,0]`
  - `E03`. project is under- or over-scoped from the beginning, so the team spends time idling or does not know where to start `[3,C,1]`, (`[13,A,3]`, `[13,B,2]`, `[13,C,4]`) -> `[13,ABC,5]`, `[15,C,1]`
- `SC2`. only analytical or documentational artifacts for a long time
  - `E04`. Development takes place only at the end of an iteration, sometimes rushed (opposite of `ASC08`) `[3,C,1]`, `[6,A,1]`, `[14,A,3]`
- `SC3`. relatively short period towards project end with a sudden increase in development efforts (i.e. rock-edge burndown, especially when viewing implementation tasks only)
  - `E05`. team rushes to deliver at least an MVP to meet the final deadline `[13,A,3]`
- `SC4`. little testing/QA and project progress tracking activities during the development period
  - `E06`. sometimes caused by improper usage of project management tools, for example logging time only at the end of a phase `[3,C,2]`
  - `E07`. too much focus on "visible" progress by managerial decision while testing/QA is neglected, which leads to an accumulation of technical debt in the long-term (Half Done Is Enough) `[12,A,1]`
- `SC5`. final product with poor code quality, many open bug reports, poor or patchy documentation
  - `E08`. too meet the final delivery date, the product quality is decreased by skipping, e.g., features or proper Q/A (alternatively, the may be delivered late, but as agreed) (`[4,B,1]`, `[4,C,2]`) -> `[4,BC,3]`, `[6,C,2]`, `[7,B,2]`, (`[13,A,4]`, `[13,C,2]`) -> `[13,AC,5]`
- `SC6`. if points `SC3` through `SC5` do not apply, (likely) the project schedule or scope is compromised (i.e., either delayed delivery or descoping occurs)
  - `E09`. team accepts change requests, re-prioritization of existing or new issues within a phase (e.g., after the start of a sprint); improper change management process `[1,A,1]`, `[3,C,4]`, (`[5,A,1]`, `[5,B,1]`) -> `[5,AB,2]`, `[9,A,1]`, `[10,B,2]`, `[12,A,3]`, (`[13,A,4]`, `[13,B,4]`) -> `[13,AB,5]`
  - `E10`. planned work is not completed and overflows into the next phase (e.g., sprint), due to, e.g., an over-challenged team (opposite of `ASC28`), misestimation, or unequal work distribution `[3,C,4]`, `[7,C,3]`, `[9,C,3]`, `[10,C,1]`, `[14,A,1]`, `[15,B,1]`
  - `E11`. iterations are too short and are artificially prolonged, forcing the team to do overtime or to truncate the workload `[3,C,4]`, `[4,C,2]`, `[7,C,1]`
- `SC7`. stark contrast between interlevel communication in project hierarchy (management - developers) during the first period (close to silence) and after realizing the problem (panic and loud noise)

Here is a list of new, empirical symptoms and causes:

- `ESC1`. poor communication (e.g., unresponsive, relayed, large overhead, or underqualified decision-maker) between stakeholders (e.g., customer) and the development team
  - `E12`. unresponsive customer or unsatisfactory (e.g., late, incomplete, or slow) communication (critical infos or materials not provided timely) (`[3,A,4]`, `[3,B,4]`) -> `[3,AB,5]`, (`[4,A,3]`, `[4,B,4]`, `[4,C,3]`) -> `[4,ABC,5]`, `[5,C,1]`, (`[6,A,3]`, `[6,B,3]`) -> `[6,AB,4]`, (`[7,A,3]`, `[7,C,3]`) -> `[7,AC,4]`, (`[10,A,3]`, `[10,B,3]`, `[10,C,2]`) -> `[10,ABC,5]`, `[11,A,1]`, `[12,A,1]`, `[13,C,1]`, `[14,A,1]`, (`[15,A,1]`, `[15,C,1]`) -> `[15,AC,2]`
  - `E13`. requirements cannot be clearly negotiated or are ambiguous `[3,C,3]`, (`[10,A,3]`, `[10,B,3]`) -> `[10,AB,4]`, (`[13,A,3]`, `[13,C,4]`) -> `[13,AC,5]`
  - `E14`. post-negotiation misunderstandings (without proper re-negotiation) `[3,A,2]`
  - `E15`. tacit misunderstanding (stakeholder and team believe they are on the same page, but they are not in actuality) `[3,A,3]`, `[4,C,1]`
  - `E16`. customer interferes with project management without properly communicating the made changes, which directly translates into a project risk `[9,B,2]`
- `ESC2`. high project risk (opposite of `ASC26`), often manifests itself through, e.g., unrealistic work item estimates, the absence of proper testing (opposite of `ASC09`), or improper documentation
  - `E17`. business requirements (tacitly/unknowingly) misinterpreted (`[3,C,3]`, `[3,A,3]`) -> `[3,AC,4]`, `[4,C,3]`, `[6,A,2]`, `[12,A,1]`
  - `E18`. finalized work does not conform to the defined specification/expectation `[3,C,2]`
  - `E19`. new functionality introduces bugs, and not enough slack was allocated during planning for the fixing (or preventing by proper testing) of these `[3,C,1]`, `[4,A,2]`
  - `E20`. tasks are done in the wrong order (Cart Before Horse), esp. development before properly analyzing and planning (`[6,A,2]`, `[6,B,2]`) -> `[6,AB,3]`, `[12,A,1]`, `[15,A,1]`
  - `E21`. imbalanced activities at the beginning, end, or during the project, such as too much focus on development early and requirements analysis later that leads to descoping, for example (`[6,B,2]`, `[6,C,2]`) -> `[6,BC,3]`, (`[7,A,1]`, `[7,B,1]`) -> `[7,AB,2]`, `[9,A,1]`, `[12,C,1]`, (`[13,A,3]`, `[13,B,2]`) -> `[13,AB,4]`
  - `E22`. lack of experience that leads to misestimation of work items `[6,C,1]`, `[9,A,1]`, `[13,C,2]`
  - `E23`. work items or goals not properly defined, absent, or defined too late `[7,C,4]`, `[13,A,4]`
  - `E24`. management fails to ascertain that the development team is available to its planned capacity (e.g., it allows the team to be affected by external factors), which has a negative impact on the progression of the project, such as descoping, quality regression, or delayed delivery `[9,C,1]`, `[11,A,1]`, `[13,A,1]`
  - `E25.` Strong dependency between stakeholders and the development team, such that the team cannot proceed very long or at all by themselves (opposite of `ASC15`) `[10,C,1]`
  - `E26`. technical difficulties in the environment, such as the infrastructure, which cause unexpected delays to, e.g., the development or deployment (`[12,A,3]`, `[12,B,1]`) -> `[12,AB,3]`, (`[15,A,1]`, `[15,C,1]`) -> `[15,AC,2]`
  - `E27`. frequent project schedule adaptions, manifested by excessive use of administrative tools, leading to a low-quality product `[13,B,4]`, `[15,A,2]`
- `ESC3`. poor usage of project management tools and methodologies which gives rise to management misinterpreting the progress and state of the project
  - `E28`. too-large goals that were not properly broken down into smaller issues `[3,C,1]`
  - `E29`. mislabeling of items; for example, marking an Epic as a Task `[4,C,1]`
  - `E30`. a discrepancy between the defined work and the actual work exists, due to, e.g., time not logged properly, issues not defined, or work completed during undocumented overtime `[3,C,1]`, `[7,C,4]`, `[10,A,1]`, `[15,A,1]`
  - `E31`. information mismanagement, for example, by duplication or using too many different systems for storing and disseminating information `[4,B,1]`

## Anti-Symptoms and -Consequences

While the absence of evidence does not mean that a Fire Drill is not present, here we gather some empirical evidence for the absence.
The symptoms and consequences below were gathered by analyzing the raters' notes in the dataset [[HON'23]](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/References.md).

- `ASC01`. communication and collaboration with the customer is seamless
- `ASC02`. no descoping which typically happens towards the end of a phase (sprint, milestone, etc.)
- `A0SC3`. timely product delivery according to agreed-upon quality
- `ASC04`. satisfaction among all stakeholders (e.g., team, customer, etc.)
- `ASC05`. regular and successful iteration evaluations that do not result in the unveiling of (large/additional) problems
- `ASC06`. clear understanding of the requirements and resulting unproblematic execution
- `ASC07`. equal (or almost equal) work distribution among team members (also: fair work distribution among differently-skilled/-tasked team members)
- `ASC08`. linear burn-down (i.e., done work is distributed uniformly, instead of, e.g., at the end of a phase)
- `ASC09`. product tested properly (e.g., appropriate tests and/or good coverage), as well regularly/continuously
- `ASC10`. start to implementing features right from the start of the project (clear requirements)
- `ASC11`. proper allocation of project resources (esp. time)
- `ASC12`. proper (planning of) distribution of time (spent) across the required activities (e.g., enough time spent on defining requirements properly)
- `ASC13`. appropriate prioritization of activities when resources (often time) become (temporarily) scarce
- `ASC14`. successful intermediate and final product deliveries (according to customer's acceptance criteria)
- `ASC15`. team can proceed at least short-term even if the customer is unavailable (good internal crisis management)
- `ASC16`. accurate work-item estimates (time, points, etc.), esp. no over-estimation (which indicates high level of uncertainty and, therefore, risk)
- `ASC17`. project management tool(s) used accordingly; e.g., proper usage of primitives (item types), the Scrum/Kanban board (or swimlanes), regular updates (all these indicate proper management)
- `ASC18`. regular activities according to used methodology (e.g., Scrum), such as daily meetings, retrospectives, and milestones
- `ASC19`. change requests, re-prioritization of existing or new issues are rejected by the team once the phase (e.g., sprint) started in which they were planned for (as should be)
- `ASC20`. proper communication among team members; direct messaging, as well as dedicated channels and often the usage of bots (from, e.g., a CI pipeline)
- `ASC21`. efficient communication with customer (that is, direct (no relays), quick, unproblematic, of high quality, tending to the necessary aspects of the product (low overhead))
- `ASC22`. stable team (no developer churn) and harmony among members
- `ASC23`. mutual understanding: effort estimations between customer and team are similar (customer understands technical challenges and team understands business requirements)
- `ASC24`. activities in right order (e.g., analysis before design before implementation etc.)
- `ASC25`. progress is reflected empirically (objectively), i.e., provably no discrepancy between reported and actual progress exists
- `ASC26`. proper risk management through, e.g., the development of prototypes
- `ASC27`. the scope may change and adapt over the course of the project (due to the agile nature), but it does not increase/widen without additional resources
- `ASC28`. team is not undersized for the project: no (steady) accumulation of non-finished work items into the next phase
- `ASC29`. when external forces and (un)forseeable events happen, development is suspended and/or the product is delayed accordingly, allowing the team to catch up (rather than forcing them to do overtime)

Sometimes, there are also signs for the __opposite__ of a Fire Drill (or other patterns), such as (non-exhaustive list):
- team underchallenged (quick completion of work items, delivery before deadline, product too simple, insufficient Q/A, etc.)
- skipping over of planning and/or analysis with direct start of implementation

In such cases, one or more other (anti-)pattern(s) perhaps apply.

## Symptoms and Consequences (in source code)

These are just some logical deductions made from known symptoms and consequences and not part of the original Fire Drill description. However, we include them here as they might still be useful.

 - rock-edge burndown of esp. implementation tasks mean there are no or just very few adaptive maintenance activities before the burning down starts
 - the long period at project start translates to few modifications made to the source code, resulting in fewer commits (lower overall relative frequency)
 - likewise, documentational artifacts have a lower _source code density_, as less functionality is delivered; this density should increase as soon as adaptive activities are registered
 - the short period at project end is characterized by a higher frequency of higher-density implementation tasks, with little to no perfective or corrective work
 - at the end of the project, code quality is comparatively lower, while complexity is probably higher, due to pressure excerted on developers in the burst phase

## Causes

- `C1`. management does not take seriously development effort (time) estimates
- `C2`. management absorbed in "various technopolitical issues (...) prevent[ing] the development staff from making progress"
- `C3`. team is happy to produce artefacts early in the project
- `C4`. requirements are complex and their prioritization is not forced early on
- `C5`. team overseeing the need to prioritize "working code over comprehensive documentation" 
- `C6`. management wants to appear the project to be on track
- `C7`. management believes it is more important to deliver complete functionality than good quality
- `C8`. project tracking and oversight is loose, easily lulled inco complacency by easy-to-reach outcomes

## (Refactored) Solution

- `RS1`. force the team to start delivering (parts of) the "consumable solution" early, possibly alongside the analysis and planning artefacts, by instituting strong project tracking and oversight related to actual outcomes
- `RS2`. it helps to follow an iterative process, architecture-driven development, and have a well-performing product owner 

## Variations (optional) 

## Example(s) (optional) 

## Related Anti-patterns

|Anti-pattern  | Relation |
|--|--|
| [Analysis Paralysis](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/catalogue/Analysis_Paralysis.md) | potential cause |
| [Collective Procrastination](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/catalogue/Collective_Procrastination.md) | more generic case |

## Notes (optional)

## Sources
[[HON'23]](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/References.md), [[SIL'15]](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/References.md), [[SOU'18]](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/References.md) [Fire Drill](https://sourcemaking.com/antipatterns/fire-drill), [[BRO'98]](https://github.com/MrShoenel/Software-process-antipatterns-catalogue/blob/7a4d8/References.md)
