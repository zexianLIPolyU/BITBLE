clc ; clear ; close all ;
% add the path of package QCLAB https://github.com/QuantumComputingLab/fable
addpath( 'QCLAB' ) ;

%% Define a matrix A in $\mathbb{C}^{2^n \times 2^n}$ and setting for BITBLE
n = 1 ;
A = randn(pow2(n),pow2(n)) + randn(pow2(n),pow2(n)) .*1j;

offset = 0 ;
logging = true ;
compr_type = 'cutoff' ;%'percentage'; 
compr_val = 1e-8 ;
circuit_sim = true ;

%% Construct the paramter of binary tree block encoding (BITBLE) using recursive multiplexed rotations
fprintf("\n\nBITBLE Block Encoding using recursive multiplexed rotations\n") ;
fprintf("------------------------------------------------------------ \n");
fprintf("parameter computing... \n") ;
tic;
[circuit1, normalization_factor, info1] = bitble( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
toc;
info1.circ

% Simulate the quantum circuit 
fprintf("matrix simulting... \n");
tic;
U_bitble1 = circuit1.matrix ;
toc;
UA_bitble1 = U_bitble1( 1:pow2(n), 1:pow2(n) ) ;
fprintf( "Frobenius norm of A = %f \n", norm(A, 'fro') ) ;
fprintf( "normalization_factor = %f \n", normalization_factor ) ;
fprintf( "norm(normalization_factor .* UA_bitble - A) = %e \n", norm(normalization_factor.*UA_bitble1 - A) ) ;

% draw quantum circuit in a Latex format
file = fopen('circuit_bitble1.tex','w');
circuit1.toTex(file);


%% Construct the paramter of binary tree block encoding (BITBLE2) using permutative multiplexed rotations
fprintf("\n\nBITBLE2 Block Encoding using permutative multiplexed rotations\n");
fprintf("------------------------------------------------------------ \n");
fprintf("parameter computing... \n");
tic;
[circuit2, normalization_factor, info2] = bitble2( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
toc;
info2.circ

% Simulate the quantum circuit 
fprintf("matrix simulting... \n");
tic;
U_bitble = circuit2.matrix ;
toc;
UA_bitble = U_bitble( 1:pow2(n), 1:pow2(n) ) ;
fprintf( "Frobenius norm of A = %f \n", norm(A,'fro') ) ;
fprintf( "normalization_factor = %f \n", normalization_factor ) ;
fprintf( "norm(normalization_factor.*UA_bitble - A) = %e \n", norm(normalization_factor.*UA_bitble - A) ) ;

% draw quantum circuit in a Latex format
file = fopen('circuit_bitble2.tex','w');
circuit2.toTex(file);


%% Construct the paramter of binary tree block encoding (BITBLE3) with normalization factor of S_q(A) using recursive multiplexed rotations
fprintf('\n\nBITBLE3 Block Encoding with normalization factor √(maxₖ‖A(k,:)‖₂ₚ₁ * maxₖ‖A(:,k)‖₂ₚ)\n');
fprintf("------------------------------------------------------------ \n");
fprintf("parameter computing... \n") ;
p = 0.5 ;
tic;
[circuit3, subnormalization_factor3, info3] = bitble3( A, p, compr_type, compr_val, logging, offset, circuit_sim ) ;
toc;
info3.circ

% Simulate the quantum circuit 
fprintf("matrix simulting... \n");
tic;
U_bitble3 = circuit3.matrix ;
toc;
UA_bitble3 = U_bitble3( 1:pow2(n), 1:pow2(n) ) ;
fprintf( "Frobenius norm of A = %f \n", norm(A,'fro') ) ;
fprintf( "normalized_factor = %f \n", subnormalization_factor3 ) ;
fprintf( "norm(normalized_factor.*UA_tbble - A) = %e \n", norm(subnormalization_factor3.*UA_bitble3 - A) ) ;

% draw quantum circuit in a Latex format
file = fopen('circuit_bitble3.tex','w');
circuit3.toTex(file);
