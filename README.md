#header
Root-flipped multiband pulses with inherently aligned echoes (ISMRM 2017 abstract nr 3955)

These are based on the Root-flipped refocusing pulses by Sharma et al (MRM 2016), modified to form spin-echoes which are inherently aligned.
This means that spin-echoes arrive (or reach their maximum) almost simultaneously, and have negligible difference in their T2 weighting.

This code was further developed from the code released by Sharma et al for the root-flipping publication, which is available at:
http://www.vuiis.vanderbilt.edu/~grissowa/

Other things that MRI physicists might find useful:
* The repo includes a script and a Bloch-simulators which is able to simulate B0-inhomogeniety (i.e. spin-dephasing due to spin-spin interactions) 
based on Adult brain white matter at 3T.
* For those working in Simultaneous Multi-Slice (SMS) Imaging, the repo also includes a function which for a multiband spin-echo (excitation + refocusing) calculates the gradient rewind area which minimizes the maximum possible in-slice phase error in radians using an fminsearch approach.


This package requires the following:
*  CVX from (http://cvxr.com/cvx/).
*  Pauly's RF tools from (http://rsl.stanford.edu/research/software.html).

11/04/2017
For contact:
Samy Abo Seada
samy.abo_seada@kcl.ac.uk

Shaihan Malik
shaihan.malik@kcl.ac.uk
King's College London, 2017.