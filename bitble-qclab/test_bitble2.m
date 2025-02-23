%clc; clear; close all;
% add the path of package QCLAB https://github.com/QuantumComputingLab/fable
addpath( 'QCLAB' );

%% Define a matrix A in $\mathbb{C}^{2^n \times 2^n}$. 
n = 3;
A = randn( pow2(n), pow2(n) ) + 1i.*randn(pow2(n),pow2(n));


%% Construct the paramter of binary tree block encoding with less CNOT gates (BITBLE2) 
offset = 0 ;
logging = 1 ;
compr_type = 'cutoff' ;%'percentage'; 
compr_val = 1e-8 ;
circuit_sim = true ;

fprintf("\n\nBITBLE2 Block Encoding \n");
fprintf("------------------------------------------------------------ \n");
fprintf("parameter computing... \n");
tic;
[circuit, normalized_factor, info] = bitble2( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
toc;
info.circ

%% Simulate the quantum circuit 
% {
fprintf("matrix simulting... \n");
tic;
U_bitble = circuit.matrix ;
toc;
UA_bitble = U_bitble( 1:pow2(n), 1:pow2(n) ) ;
fprintf( "Frobenius norm of A = %f \n", norm(A,'fro') ) ;
fprintf( "normalized_factor = %f \n", normalized_factor ) ;
fprintf( "norm(normalized_factor.*UA_tbble - A) = %e \n", norm(normalized_factor.*UA_bitble - A) ) ;
%}

