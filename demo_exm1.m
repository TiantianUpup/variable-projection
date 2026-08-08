%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo_exm1 reproduces the numerical experiments in Subsection 6.1
% (Tables 1-3).
%
% This example considers tensors with a large third dimension p.
%
% The experiments reported in the paper consider
% p = 100, 1000, and 10000 with different values of m and n.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear all; close all;

addpath('examples');

% Case 1: corresponds to Table 1
% m=100;
% n=100;
% p=100;

% Case 2: corresponds to Table 2
m=30;
n=40;
p=1000;

% Case 3: corresponds to Table 3
% m=100;
% n=100;
% p=10000;

% approximation rank
r=10;
test_exm1(m,n,p,r)

% case 2
% m=30;
% n=40;
% p=1000;

% case 3
% m=100;
% n=100;
% p=10000;