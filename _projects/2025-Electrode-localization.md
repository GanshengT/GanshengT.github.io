---
title: "Electrode-localization"
collection: projects
type: "innovation project"
permalink: /projects/2025-Electrode-localization
venue: "Department of Biomedical Engineering"
date: 2025-01-28
location: "St. Louis, USA"
---




Algorithms [![DOI](https://zenodo.org/badge/891756132.svg)](https://doi.org/10.5281/zenodo.14217838) in MATLAB for **objectively and automatically** localizing intracranial contac, based on CT images, MRI images, and freesurfer segmentation.

Workflow
============================

Implementation
--------

|  | Entry                                                        | Description                                                 |
|--| --------                                                     |------------------------------------------------------------ |
|01| Initialization | Use predictive algorithms or manually select initial electrode locations based on CT      |
|02| Searching | Use graph-based search with brightness threshold going from high to low until a connected volume is identified  |
|03| Validation | The center of the identified bright cloud is used as electrode centroid. These locations will be exported to freeview for validation |
|04| Customization | If the electrode are evenly spaced, based on shank direction, we can initialize without manual selection|




Future direction
---------
Achieve non-parametric optimization to maximize dependence on subjective choice.


Challenges 
---------


