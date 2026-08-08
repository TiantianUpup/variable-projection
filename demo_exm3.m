%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo_exm3 reproduces the numerical experiments in Subsection 6.3
% (Table 5 and Table 6).
%
% The factor matrix W is ill-conditioned.
% The parameter kappa controls the condition number of W.
%
% The experiments reported in the paper use kappa = 1e-6 and 1e-8.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear all; close all;
addpath('examples');

m=30;
n=40;
p=50;

% approximation rank
r=10;

% Case 1: corresponds to Table 5
% kappa=1e6;

% Case 2: corresponds to Table 6
kappa=1e6;
test_exm3(m,n,p,r,kappa);

% case 2
% m=30;
% n=40;
% p=1000;

% case 3
% m=100;
% n=100;
% p=10000;