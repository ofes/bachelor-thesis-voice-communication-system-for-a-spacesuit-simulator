%% main: Genetic algorithm 
% In this short program the use of the genetic algorithm is being tested.

clc;
clear all;
close all;

FitFcn = @myFitness;        % Giving the fitnes function another name.
nvars = 2;                  % Number of variables of the fitness function.

[x, fval] = ga(FitFcn, nvars);

x
fval