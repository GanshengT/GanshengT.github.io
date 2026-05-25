---
title: "Explainable machine learning for cognitive neuroscience"
collection: projects
type: "research framework"
permalink: /projects/2025-explainable-machine-learning-cognitive-neuroscience
venue: "Washington University in St. Louis"
date: 2025-11-15
location: "St. Louis, USA"
summary: "An explainable machine learning framework for linking neural dynamics to cognitive behavior while preserving feature-level scientific interpretation. The project uses interpretable models to identify neural patterns that support visual encoding and candidate markers for closed-loop cognitive neurotechnology."
related_research_label: "Brain-body coupling"
related_research_url: "/research/#brain-body-coupling"
related_publications:
  - label: "Mind's eye: Saccade-related evoked potentials support visual encoding in humans"
    url: "/publication/2025-11-15-minds-eye-saccade-related-evoked-potentials-support-visual"
---

This project develops an explainable machine learning framework for cognitive neuroscience. The goal is to learn complex relationships between neural signals and cognitive behavior, then explain model decisions in a way that supports scientific interpretation rather than only predictive performance.

Framework
============================

The current workflow uses random forests to model nonlinear relationships between neural features and behavioral outcomes. We then use Shapley-value-based explanations to quantify how each feature contributes to model predictions at the individual and population levels.

Interpretability
--------

The explanation layer follows the Shapley-value formulation described in [Interpretable Machine Learning by Christoph Molnar](https://christophm.github.io/interpretable-ml-book/shapley.html). In this project, that framework is used to identify which neural features are most informative for prediction and how those contributions vary across behavioral states or task conditions.

Current application
--------

One target application is visual encoding during natural eye movements. In the saccade-related evoked potential work, the framework is intended to connect multichannel neural dynamics to cognitive behavior while preserving an interpretable account of feature importance.

Next steps
--------

- Expand the feature space to include temporal, spectral, and cross-channel interaction terms.
- Compare random forests with other nonlinear models while keeping the explanation layer scientifically interpretable.
- Use Shapley profiles to identify candidate neural markers for closed-loop cognitive neurotechnology.
