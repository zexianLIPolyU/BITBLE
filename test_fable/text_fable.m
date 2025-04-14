clc; clear; close all;
% add the path of package QCLAB
addpath( 'QCLAB' );

%% Define a matrix A in $\mathbb{C}^{2^n \times 2^n}$. 
n = 4;
A = randn( pow2(n), pow2(n) ) + 1i.*randn(pow2(n),pow2(n));


%% Construct the paramter of binary tree block encoding (BITBLE) 
offset = 0 ;
logging = 1 ;
compr_type = 'cutoff' ;%'percentage'; 
compr_val = 1e-8 ;
circuit_sim = true ;

fprintf("\n\nFABLE Block Encoding \n");
fprintf("------------------------------------------------------------ \n");
fprintf("parameter computing... \n");
tic;
[circuit, OA, alpha, info] = fable( A, compr_type, compr_val, logging, circuit_sim) ;
toc;

%% Simulate the quantum circuit 
fprintf("matrix simulting... \n");
% tic;
% U_bitble= circuit.matrix ;
% toc;
% UA_bitble = U_bitble( 1:pow2(n), 1:pow2(n) ) ;
% normalized_factor = real(A(1,1)) ./ real(UA_bitble(1,1)) ;
% fprintf( "norm(normalized_factor.*UA_tbble - A) = %e \n", norm(normalized_factor.*UA_bitble - A) ) ;
subnormalized_factor = pow2(n) * max(max(abs(A))) ;
