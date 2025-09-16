# MLM - Dynamic Survival Modelling of Risk Credit Default
This study was undertaken in the summer of 2025 as part of the Mary Lister McCammon Summer Fellowship, under the supervision of F. Javier Rubio from the Department of Statistical Science at UCL.

This repository contains the code and documentation for the analysis. The objective of this work is to examine the dynamics of mortgage loans in the period preceding default. To this end, we adapt survival models originally developed in biostatistics to the context of the mortgage market. Within this framework, the event is defined as the granting of a mortgage loan, whereas death is reinterpreted as the point at which the loan defaults. This conceptual shift provides a novel and pioneering contribution to the existing literature.


The repository is organized as follows:

1. Simulation.R – Generates simulated data for the UK mortgage sector.

2. StaticModelling.R – Implements Proportional Hazards and Accelerated Failure Time (AFT) models using the simulated data.

2. DynamicModelling.R – Implements the Landmarking model using the simulated data.

3. Prediction.R – Compares predicted probabilities at the individual level.

4. Testing.R – Provides a pseudo-evaluation metric for the AFT and Landmarking models.


A presentation summarizing this work can be accessed here: https://github.com/anamariapatron/MLM-On-The-Applications-Of-Survival-Modelling-In-The-Mortgage Sector/blob/fba55ceab0114b038a0114d465eefc3a30d17bac/documents/Presentation_MLM__Survival_Modelling_applied_to_mortgage_sector.pdf

A poster of this work can be found here: https://github.com/anamariapatron/MLM-On-The-Applications-Of-Survival-Modelling-In-The-Mortgage-Sector/blob/fba55ceab0114b038a0114d465eefc3a30d17bac/documents/Poster_AMPP_VF.pdf

Mary Lister McCammon Summer Fellowship is a funded opportunity developed at Imperial College London and UCL. More information here: https://www.imperial.ac.uk/mathematics/postgraduate/the-mary-lister-mccammon-summer-research-fellowship/
