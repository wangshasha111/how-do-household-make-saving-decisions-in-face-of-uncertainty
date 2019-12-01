%% Research Project 
% Economics 712
% Fall 2019
% Dirk Krueger

% Value function iteration taken from "Basic RBC model with full depreciation" by Jesus Fernandez-Villaverde Haverford, July 3, 2013

% 2019-11-6 22:06:22


clear;
close all;

cd 'E:\Dropbox\fall 19-20\dirk_krueger\project'

%% Part I Partial Equilibrium

% exercise 2
ssigma = 1;
ddelta = 0.8;
ssigmaY = [0.2,0.4];
rho = 0.04;
r = 0.02;
M = 1000; % number of simulated paths
T = 61; % 61-period simulated paths

