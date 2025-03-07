clc;clear;close all
% csd()  cosine-sine decomposition 

addpath("iwoodsawyer-csd-a23bac9"); % loading csd() % https://www.mathworks.com/matlabcentral/fileexchange/50402-cosine-sine-decomposition
addpath("QCLAB");  %https://github.com/QuantumComputingLab/fable

%% Define a matrix A in $\mathbb{C}^{2^n \times 2^n}$ and setting for SIABLE
n = 3 ;
m = pow2(n) ;
A = randn(m,m) +randn(m,m).*1i ;

offset = 3 ;
logging = 1 ;
compr_type = 'cutoff' ;%'percentage'; 
compr_val = 1e-8 ;
circuit_sim = true ;

% Simulate the quantum circuit 
fprintf("\nSIABLE Block Encoding");
fprintf("\n------------------------------------------------------------ \n");
fprintf("parameter computing... \n") ;
tic ;
[circuit, subnormalized_factor, info] = siable( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
toc ;
fprintf("1.0001 * 2-norm of A = %f \n",1.0001 *norm(A,2)) ;
fprintf("normalized_factor = %f \n",subnormalized_factor) ;
M1 = circuit.matrix;
fprintf("norm(subnormalized_factor.*M1(1:m,1:m)-A) = %e \n",norm(subnormalized_factor.*M1(1:m,1:m)-A)) ;
if logging, info.circ; end 


