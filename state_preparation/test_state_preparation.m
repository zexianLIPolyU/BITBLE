clc; clear; close all;
addpath('QCLAB');

logging = false; % no record
circuit_sim = true ;
N = 8;
n = log2( N ) ;


%% State Praprartion with real amplitude

state_real = randn(N,1) ;
state_real = state_real ./ norm(state_real) ;
[circuit1, info1] = binary_tree_statepreparation(state_real, logging, circuit_sim) ;
% output
disp("Real State Preparation:");
disp("circuit:");
circuit1.draw() ;
disp("original state:");
disp(state_real);
mat = circuit1.matrix ;
disp("state preparation:");
disp(mat(:,1));
fprintf("norm(original state - prepared state) = %e\n",norm(state_real - mat(:,1)));
disp("-----------------------------------------------------------------------");


%% State Praprartion with complex amplitude
state_complex = randn(N,1) + randn(N,1).*1i ;
state_complex = state_complex ./ norm(state_complex) ;
[circuit2, info2] = binary_tree_statepreparation(state_complex, logging, circuit_sim) ;
disp("Real State Preparation:");
disp("circuit:");
circuit2.draw() ;
disp("original state:");
disp(state_complex);
mat2 = circuit2.matrix ;
disp("state preparation:");
disp(mat2(:,1));
fprintf("norm(original state - prepared state) = %e\n",norm(state_complex - mat2(:,1)));
disp("-----------------------------------------------------------------------");

