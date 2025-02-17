# ScBM-VAR
ScBM-VAR is a dynamic network model in which the network structure in each specific period, called the season, follows a degree-corrected Stochastic co-BlockModel (ScBM), and the observations are collected by variants of Vector AutoRegressive models (VAR). The Periodic Vector AutoRegressive (PVAR) and Vector Heterogenous AutoRegressive (VHAR) models are considered in the study, and the estimated transition matrices form the initial network structures of each season. Once the transition matrices are given, both left and right singular vectors of all matrices are smoothed by PErsistent Communities by Eigenvector Smoothing (PisCES) algorithm, and the normalized singular vectors are used for spectral clustering. A fully functioning package is under development and is expected to be released soon. The corresponding manuscript is Kim and Baek (2024).
