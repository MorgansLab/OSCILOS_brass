---
title: 'OSCILOS_brass: an acoustic solver for brass instruments'
bibliography: paper.bib
tags:
- brass instruments
- acoustics
- aeroacoustics
- musical acoustics
- MATLAB
authors:
- name: Alexander MacLaren
  orcid: 0000-0002-5835-216X
  affiliation: 1
- name: Aimee S. Morgans
  affiliation: 1
- name: Renaud Gaudron^[Corresponding author, r.gaudron\@imperial.ac.uk]
  affiliation: 1
affiliations:
- name: Department of Mechanical Engineering, Imperial College London, UK
  index: 1
date: 13 March 2021
---

# Summary

 The sound, or timbre, of musical instruments is governed by the 
 distribution of their modal frequencies. In brass instruments like the 
 trumpet, french horn and trombone, the geometry of the internal bore 
 is the principal factor determining the modal frequencies. Harmonic 
 analysis of brass instruments is necessary to assist bore shape 
 optimisation by instrument makers, and to enable research into lip 
 dynamics, the origins of 'brassiness', and other prevailing questions 
 in brass wind acoustics. OSCILOS_brass provides an easy-to-use 
 platform for determining the modal frequencies of a given instrument 
 geometry, be it trombone or drainpipe.


# Statement of need

OSCILOS_brass is intended as a tool to assist brass wind acoustics 
research, allowing the effects of nuances in instrument geometry on 
sound and playability to be explored. OSCILOS_brass is based on 
OSCILOS_long [@Li2017], an open source code for simulating 
combustion instabilities, whose predictions have been validated against experiments [@Han2015; @Xia2019]. OSCILOS_brass is written in MATLAB, is highly modular, and 
is straightforward to run and edit. It represents an instrument bore as 
a sequence of connected cylindrical finite elements, using a 1-D plane 
wave approximation. A variety of inlet and exit acoustic boundary 
conditions are available, including open, closed, Levine-Schwinger and 
user-defined boundary conditions. The mean flow is calculated assuming 
1-D flow conditions, with changes only across element interfaces. 
Extensive documentation is provided in the accompanying technical 
report and user guide [@MacLaren2021].

A typical tenor trombone geometry is shown in \autoref{fig:geom}.

![Trombone bore profile from [@Bilbao2013] as represented by 
OSCILOS_brass.\label{fig:geom}](figures/TromboneGeometry.png)

OSCILOS_brass calculates Equivalent Fundamental Pitch (EFP) deviation 
defined by \autoref{eqn:EFP} [@Chick2004], where 
$f_i$ is the frequency of the $i$th mode, and $F$ is the reference 
fundamental pitch, conventionally taken as $f_4/4$, the note to which 
brass players often tune their instruments. The unit of this definition 
is the \verb|[cent]|, equal to $1/100$th of a semitone, or a frequency 
ratio of $\sqrt[1200]{2}$.

\begin{equation}\label{eqn:EFP}
	\mathrm{EFP}(f_i) = \frac{1200}{\log(2)}\log\left[\frac{f_i}{iF}\right]
\end{equation}

Good agreement with published measurements and simulations for this 
geometry, and for a tuba geometry from a separate study [@Norman2013], 
is demonstrated in \autoref{fig:Tromb}.

![OSCILOS_brass EFP comparison to results for trombone geometry from 
[@Bilbao2013] (left), and to results for tuba geometry from 
[@Norman2013] (right).\label{fig:Tromb}](figures/TromboneTubaEFP.png)

EFP is shown for 4 combinations of boundary 
conditions in \autoref{fig:trombBC}, as calculated for 
the trombone geometry in \autoref{fig:geom} by OSCILOS_brass.

![EFP output by OSCILOS_brass for 4 \[Inlet - Outlet\] sets of boundary 
conditions applied to the trombone geometry from 
[@Bilbao2013] (left), and to results for tuba geometry from 
[@Norman2013] (right) \label{fig:trombBC}](figures/TromboneTubaBCsEFP.png)

# Acknowledgements

This work was supported in part by the European Research Council (ERC) 
Consolidator Grant AFIRMATIVE (2018–2023). A. MacLaren acknowledges the 
Mechanical Engineering UROP Bursary 2020 from Imperial College London.

# References
