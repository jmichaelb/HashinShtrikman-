# HashinShtrikman-
Matlab Code associated with "Determination of Hashin-Shtrikman Bounds on the Isotropic Effective Elastic Moduli of Polycrystals of any Symmetry"
Two files are included in this repository:

test_HSBounds.m
HSBounds.m

The first contains example data, makes the call to the function, and includes a listing of expected results
The second is the implementation determining the Hashin-Shtrikman Bounds:

function [hs,vrh,ko_go,ustart]=HSBounds(cij)
% This function finds the Hashin-Shtrikman bounds for material with elastic
% moduli cij in a 6x6 matrix.  It is the implementation of the variational 
% equations first derived by Hashin and Shtrikman in 1963 and later
% articulated in papers by Peselnick and Meister and Watt and Peselnick. 
% Usage:
%       [hs,vrh,ko_go,ustart]=HSBounds(cij)
% where:
%     cij is a 6x6 matrix of elastic constants of any symmetry  (GPa units)
%     hs is a 2x2 matrix with upper and lower "optimal" bounds for K and G
%     vrh gives, for reference, the Voigt-Reuss-Hill bounds
%     ko_go gives the properties of the reference material at the optimal point
%     ustart gives the starting point used for "go" in the calculations
%
% The following "nested" custom functions are included:
%    1. [xmax,hs]=edgeu('pos' or 'neg'); finds the positive definite
%        boundary at small ko and the negative definite boundary at large ko
%    2. [k,hsl]=edgek(y,'pos' or 'neg'); finds the positive/negative
%        boundary at fixed uo with increasing ko
%    3.  y=lowerbound(x) returns the hs value at the boundaries - used in
%         search for optimal values
%    4.  y=upperbound(x) returns the hs value at the boundaries - used in
%         search for optimal values
%    5.  [hs,value]=hscalc(ko,go,cij) does the Hashin-Shrtikman
%         calculations
%    6.  [K,G]=VRH(C) calculates Voigt-Reuss-Hill bounds
% The built-in to MATLAB search function fminbnd is used to find the optimal bounds
%  JMB 2013
