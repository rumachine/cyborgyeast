Cyborg Yeast
===========

In the Nature Biotechnology publication 'In silico feedback for in vivo regulation of a gene expression circuit' (<a href="http://www.nature.com/nbt/journal/v29/n12/full/nbt.2018.html" target="_blank">here</a>, <a href="http://www.bbc.com/news/science-environment-15598887" target="_blank">here</a>, and <a href="http://control.ee.ethz.ch/index.cgi?page=publications&action=details&id=3903" target="_blank">here</a>), we showed that difficulties in regulating cellular behavior with synthetic biological circuits may be circumvented using in silico feedback control (i.e. an algorithm run on a computer). By tracking a circuit's output in Saccharomyces cerevisiae in real time, we precisely controlled its behavior using an in silico feedback algorithm to compute regulatory inputs implemented through a genetically encoded light-responsive module.

This repository contains all software used to create the mathematical model of the genetic process (i.e. the evolution of average mRNA and protein within the cell culture as a function of the light sequence applied to the cell culture), simulate the evolution of the genetic process over multiple hours, and apply feedback control algorithms to the genetic process in order to steer the average protein abundance to a desired set point.

The details of the mathematical model of the genetic process, the paramter fitting technique, and the real time closed-loop implementation are outlined in the <a href="http://www.nature.com/nbt/journal/v29/n12/extref/nbt.2018_S1.pdf" target="_blank">supplementary text</a>.

Parameter Fitting
==================

To estimate the free model parameters from data, we employed an Approximate Bayesian Computation (ABC) method. Although the maximum likelihood parameter estimates were eventually used (reported in Supplementary Table 1), the posterior distributions obtained helped us quantify the uncertainty regarding each parameter, as well as identify correlations among parameters.

The software and data used to implement the ABC method is in the Parameter Identification folder. ABCinit\_light loads the data from the experiment and initialises the optimisation structure. ABCgo\_light then runs the ABC method; the output is a distribution over the model parameter space.

Simulation
==========

A simulation of the genetic process (and corresponding trajectory data from the actual experiment) is in the Simulation (Open Loop) Experiments subfolders (for different experiments, i.e. experiments with different light sequences due to different set point objectives). The simulation can be run using CTL\_insilico.

Real-time Implementation
========================

The software used to run the feedback control experiments can be found in the Feedback (Closed Loop) Experiments subfolders (for different experiments). The file CTL implements the estimation and control algorithms on the actual data that was extracted from the cell culture in real time.

