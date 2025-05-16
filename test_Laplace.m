% testing tbble
%% L. Heisenberg Hamiltonians

clc;clear; close all
addpath( "Examples");
addpath( "QCLLB" ) ;
addpath( "bitble-qclab" ) ;
% load("record_Laplace.mat") ;

offset = 0 ;
logging = 1 ;
compr_type = 'cutoff' ;%'percentage'; 
compr_val = 1e-8 ;
circuit_sim = false ;
p = 0.5 ;




%% loading unitary for Qiskit

Laplace_mat = struct() ;
% 1D Non-periodic
for i = 1:6
    numqubit = i+1;
    dimesion = 1;
    boundary_condition = false;
    L = LaplaceOperator(numqubit, dimesion, boundary_condition);
    [BlockEncodingMatrix, normalization_factor] = BlockEncoding(L) ;
    Laplace_mat(1,i).M = BlockEncodingMatrix;
    Laplace_mat(1,i).normalization_factor = normalization_factor;
end
% 1D Periodic
for i = 1:6
    numqubit = i+1;
    dimesion = 1;
    boundary_condition = true;
    L = LaplaceOperator(numqubit, dimesion, boundary_condition);
    [BlockEncodingMatrix, normalization_factor] = BlockEncoding(L) ;
    Laplace_mat(2,i).M = BlockEncodingMatrix;
    Laplace_mat(2,i).normalization_factor = normalization_factor;
end

% 2D Non-periodic
qubits_case = [1,1; 1,2; 2,1; 2,2; 3,1; 1,3; 2,3; 3,2; 4,1; 1,4; 3,3; 3,4; 4,3] ;
dimesion = 2 ;
boundary_condition = [false,false] ;

for i = 1 : size( qubits_case,1 )
    numqubit = qubits_case(i,:) ;
    L = LaplaceOperator(numqubit, dimesion, boundary_condition) ;
    [BlockEncodingMatrix, normalization_factor] = BlockEncoding(L) ;
    Laplace_mat(3,i).M = BlockEncodingMatrix;
    Laplace_mat(3,i).normalization_factor = normalization_factor;
end

% 2D Periodic
qubits_case = [1,1; 1,2; 2,1; 2,2; 3,1; 1,3; 2,3; 3,2; 4,1; 1,4; 3,3; 3,4; 4,3] ;
dimesion = 2 ;
boundary_condition = [ true, true ] ;

for i = 1 : size( qubits_case,1 )
    numqubit = qubits_case(i,:) ;
    L = LaplaceOperator( numqubit, dimesion, boundary_condition ) ;
    [BlockEncodingMatrix, normalization_factor] = BlockEncoding(L) ;
    Laplace_mat(4,i).M = BlockEncodingMatrix;
    Laplace_mat(4,i).normalization_factor = normalization_factor;
end

save("Laplace_mat.mat","Laplace_mat") ;



% BlockEncoding()
function [BlockEncodingMatrix, normalized_factor] = BlockEncoding(A)
% Generate the (1.0001*norm(A),1,0)-block encoding 'BlockEncodingA' of 'A'
% input: A \in \mathbb{C}^{n\times n}
% output: normlized_factor = 1.0001*norm(A); A = A ./ normlized_factor;
%         BlockEncodingMatrix = [A, sqrt(eye(n)-A*A'); sqrt(eye(n)-A'*A), A']
    
    assert( size(A,1) == size(A,2) );
    n = size(A,1);
    normalized_factor = 1.0001*norm(A);
    A = A ./ normalized_factor;
    [uA12, sA12] = eig( eye(n)-A*A' );
    UA12 = uA12 * diag( sqrt(diag(sA12)) ) * uA12';
    [uA21, sA21] = svd( eye(n)-A'*A );
    UA21 = uA21 * diag( sqrt(diag(sA21)) ) * uA21';
    BlockEncodingMatrix = [ A,  UA12; -UA21, A' ];
    %  or BlockEncodingMatrix = [A,  UA12'; UA21, -A'];
end % end of BlockEncoding




record_Laplace = struct() ;
%% C. Elliptic partial differential equations - 1D Non-periodic

index = 1;
for i = 1:6
    numqubit = i+1;
    dimesion = 1;
    boundary_condition = false;
    L = LaplaceOperator(numqubit, dimesion, boundary_condition);




    fprintf("\nFor %d qubit non-periodic 1D BC", i);
    fprintf("\nHeisenberg Hamiltonians\n");
    fprintf("------------------------------------------------------------ \n");

    n = log(size(L,1))/log(2);

    fprintf("\n\nFLBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, ~, alphaL, info] = fable( L, compr_type, compr_val, logging, circuit_sim );
    time3 = toc;
    record_Laplace(index,i).fable_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).fable_nRZ = info.circ.nRZ ;
    end
    normalization_factor = max(max(abs(L))) .* pow2(n) ;
    record_Laplace(index,i).fable_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).fable_normalization_factor = normalization_factor ;
    record_Laplace(index,i).fable_time = time3 ;
    if ~real(L)
        record_Laplace(index ,i).fable_nG_metric = normalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    else
        record_Laplace(index,i).fable_nG_metric = normalization_factor .* info.circ.nRY ;
    end
    record_Laplace(index,i).fable_nCNOT_metric = normalization_factor .* info.circ.nCNOT ;

    fprintf("\n\nBITBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble( L, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time2 = toc;
    record_Laplace(index,i).bitble_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble_subnormalized_factor = subnormalization_factor ;
    record_Laplace(index,i).bitble_time = time2 ;
    record_Laplace(index,i).bitble_nG_metric = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_Laplace(index,i).bitble_nCNOT_metric = subnormalization_factor .* info.circ.nCNOT ;

    fprintf("\n\nBITBLE2 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble2( L, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time2 = toc;
    record_Laplace(index,i).bitble2_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble2_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble2_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble2_subnormalized_factor = subnormalization_factor ;
    record_Laplace(index,i).bitble2_time = time2 ;
    record_Laplace(index,i).bitble2_nG_metric = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_Laplace(index,i).bitble2_nCNOT_metric = subnormalization_factor .* info.circ.nCNOT ;

    
    fprintf("\n\nBITBLE3 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, normalization_factor, info] = bitble3( L, p, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time1 = toc;
    fprintf( "Frobenius norm of F = %f \n", norm(L,'fro') ) ;
    fprintf( "normalized_factor = %f \n", normalization_factor ) ;
    record_Laplace(index,i).bitble3_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble3_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble3_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble3_normalized_factor = normalization_factor ;
    record_Laplace(index,i).bitble3_time = time1 ;
    record_Laplace(index,i).bitble3_nG_metric = normalization_factor .* (info.circ.nRY + ~isreal(L).*info.circ.nRZ) ;
    record_Laplace(index,i).bitble3_nCNOT_metric = normalization_factor .* info.circ.nCNOT ;

save("record_Laplace.mat","record_Laplace") ;
end






%% 1D Non-periodic

index = 2;


for i = 1:6
    numqubit = i+1;
    dimesion = 1;
    boundary_condition = true;

    L = LaplaceOperator(numqubit, dimesion, boundary_condition);
    n = log(size(L,1))/log(2);
    
    fprintf("\n\nFLBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, ~, alphaL, info] = fable( L, compr_type, compr_val, logging, circuit_sim );
    time3 = toc;
    record_Laplace(index,i).fable_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).fable_nRZ = info.circ.nRZ ;
    end
    normalization_factor = max(max(abs(L))) .* pow2(n) ;
    record_Laplace(index,i).fable_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).fable_normalization_factor = normalization_factor ;
    record_Laplace(index,i).fable_time = time3 ;
    if ~real(L)
        record_Laplace(index ,i).fable_nG_metric = normalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    else
        record_Laplace(index,i).fable_nG_metric = normalization_factor .* info.circ.nRY ;
    end
    record_Laplace(index,i).fable_nCNOT_metric = normalization_factor .* info.circ.nCNOT ;

    
    fprintf("\n\nBITBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble( L, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time2 = toc;
    record_Laplace(index,i).bitble_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble_subnormalized_factor = subnormalization_factor ;
    record_Laplace(index,i).bitble_time = time2 ;
    record_Laplace(index,i).bitble_nG_metric = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_Laplace(index,i).bitble_nCNOT_metric = subnormalization_factor .* info.circ.nCNOT ;


    fprintf("\n\nBITBLE2 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble2( L, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time2 = toc;
    record_Laplace(index,i).bitble2_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble2_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble2_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble2_subnormalized_factor = subnormalization_factor ;
    record_Laplace(index,i).bitble2_time = time2 ;
    record_Laplace(index,i).bitble2_nG_metric = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_Laplace(index,i).bitble2_nCNOT_metric = subnormalization_factor .* info.circ.nCNOT ;

    
    fprintf("\n\nBITBLE3 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, normalization_factor, info] = bitble3( L, p, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time1 = toc;
    fprintf( "Frobenius norm of F = %f \n", norm(L,'fro') ) ;
    fprintf( "normalized_factor = %f \n", normalization_factor ) ;
    record_Laplace(index,i).bitble3_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble3_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble3_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble3_normalized_factor = normalization_factor ;
    record_Laplace(index,i).bitble3_time = time1 ;
    record_Laplace(index,i).bitble3_nG_metric = normalization_factor .* (info.circ.nRY + ~isreal(L).*info.circ.nRZ) ;
    record_Laplace(index,i).bitble3_nCNOT_metric = normalization_factor .* info.circ.nCNOT ;

save("record_Laplace.mat","record_Laplace") ;
end



%% 2D Non-periodic

index = 3;
qubits_case = [1,1; 1,2; 2,1; 2,2; 3,1; 1,3; 2,3; 3,2; 4,1; 1,4; 3,3; 3,4; 4,3] ;
dimesion = 2 ;
boundary_condition = [false,false] ;

for i = 1 : size( qubits_case,1 )
    numqubit = qubits_case(i,:) ;
    L = LaplaceOperator(numqubit, dimesion, boundary_condition) ;
    fprintf("\nFor (%d,%d) non-periodic 2D BC", numqubit(1), numqubit(2) ) ;

    n = log(size(L,1))/log(2);
    
    fprintf("\n\nFLBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, ~, alphaL, info] = fable( L, compr_type, compr_val, logging, circuit_sim );
    time3 = toc;
    record_Laplace(index,i).fable_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).fable_nRZ = info.circ.nRZ ;
    end
    normalization_factor = max(max(abs(L))) .* pow2(n) ;
    record_Laplace(index,i).fable_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).fable_normalization_factor = normalization_factor ;
    record_Laplace(index,i).fable_time = time3 ;
    if ~real(L)
        record_Laplace(index ,i).fable_nG_metric = normalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    else
        record_Laplace(index,i).fable_nG_metric = normalization_factor .* info.circ.nRY ;
    end
    record_Laplace(index,i).fable_nCNOT_metric = normalization_factor .* info.circ.nCNOT ;


    fprintf("\n\nBITBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble( L, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time2 = toc;
    record_Laplace(index,i).bitble_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble_subnormalized_factor = subnormalization_factor ;
    record_Laplace(index,i).bitble_time = time2 ;
    record_Laplace(index,i).bitble_nG_metric = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_Laplace(index,i).bitble_nCNOT_metric = subnormalization_factor .* info.circ.nCNOT ;
    


    fprintf("\n\nBITBLE2 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble2( L, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time2 = toc;
    record_Laplace(index,i).bitble2_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble2_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble2_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble2_subnormalized_factor = subnormalization_factor ;
    record_Laplace(index,i).bitble2_time = time2 ;
    record_Laplace(index,i).bitble2_nG_metric = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_Laplace(index,i).bitble2_nCNOT_metric = subnormalization_factor .* info.circ.nCNOT ;

    
    fprintf("\n\nBITBLE3 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, normalization_factor, info] = bitble3( L, p, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time1 = toc;
    fprintf( "Frobenius norm of F = %f \n", norm(L,'fro') ) ;
    fprintf( "normalized_factor = %f \n", normalization_factor ) ;
    record_Laplace(index,i).bitble3_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble3_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble3_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble3_normalized_factor = normalization_factor ;
    record_Laplace(index,i).bitble3_time = time1 ;
    record_Laplace(index,i).bitble3_nG_metric = normalization_factor .* (info.circ.nRY + ~isreal(L).*info.circ.nRZ) ;
    record_Laplace(index,i).bitble3_nCNOT_metric = normalization_factor .* info.circ.nCNOT ;

save("record_Laplace.mat","record_Laplace") ;
end


%% 2D Periodic
index = 4;
qubits_case = [1,1; 1,2; 2,1; 2,2; 3,1; 1,3; 2,3; 3,2; 4,1; 1,4; 3,3; 3,4; 4,3] ;
dimesion = 2 ;
boundary_condition = [ true, true ] ;

for i = 1 : size( qubits_case,1 )
    numqubit = qubits_case(i,:) ;
    L = LaplaceOperator( numqubit, dimesion, boundary_condition ) ;
    fprintf("\nFor (%d,%d) periodic 2D BC", numqubit(1),numqubit(2)) ;
    n = log(size(L,1))/log(2);
    
    fprintf("\n\nFLBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, ~, alphaL, info] = fable( L, compr_type, compr_val, logging, circuit_sim );
    time3 = toc;
    record_Laplace(index,i).fable_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).fable_nRZ = info.circ.nRZ ;
    end
    normalization_factor = max(max(abs(L))) .* pow2(n) ;
    record_Laplace(index,i).fable_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).fable_normalization_factor = normalization_factor ;
    record_Laplace(index,i).fable_time = time3 ;
    if ~real(L)
        record_Laplace(index ,i).fable_nG_metric = normalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    else
        record_Laplace(index,i).fable_nG_metric = normalization_factor .* info.circ.nRY ;
    end
    record_Laplace(index,i).fable_nCNOT_metric = normalization_factor .* info.circ.nCNOT ;


    fprintf("\n\nBITBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble( L, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time2 = toc;
    record_Laplace(index,i).bitble_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble_subnormalized_factor = subnormalization_factor ;
    record_Laplace(index,i).bitble_time = time2 ;
    record_Laplace(index,i).bitble_nG_metric = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_Laplace(index,i).bitble_nCNOT_metric = subnormalization_factor .* info.circ.nCNOT ;



    fprintf("\n\nBITBLE2 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble2( L, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time2 = toc;
    record_Laplace(index,i).bitble2_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble2_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble2_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble2_subnormalized_factor = subnormalization_factor ;
    record_Laplace(index,i).bitble2_time = time2 ;
    record_Laplace(index,i).bitble2_nG_metric = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_Laplace(index,i).bitble2_nCNOT_metric = subnormalization_factor .* info.circ.nCNOT ;

    
    fprintf("\n\nBITBLE3 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, normalization_factor, info] = bitble3( L, p, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time1 = toc;
    fprintf( "Frobenius norm of F = %f \n", norm(L,'fro') ) ;
    fprintf( "normalized_factor = %f \n", normalization_factor ) ;
    record_Laplace(index,i).bitble3_nRY = info.circ.nRY ;
    if ~isreal(L)
    record_Laplace(index,i).bitble3_nRZ = info.circ.nRZ ;
    end
    record_Laplace(index,i).bitble3_nCNOT = info.circ.nCNOT ;
    record_Laplace(index,i).bitble3_normalized_factor = normalization_factor ;
    record_Laplace(index,i).bitble3_time = time1 ;
    record_Laplace(index,i).bitble3_nG_metric = normalization_factor .* (info.circ.nRY + ~isreal(L).*info.circ.nRZ) ;
    record_Laplace(index,i).bitble3_nCNOT_metric = normalization_factor .* info.circ.nCNOT ;

save("record_Laplace.mat","record_Laplace") ;
end







%% Plot 1D Non-periodic


clc;clear; close all;
load("record_Laplace.mat") ;
load("Laplace_mat.mat") ;

MarkerSize1 = 8 ;
MarkerSize2 = 6 ;
LineWidth = 0.5 ;

index = 4;
% index = 1 -- 1D non-periodic; index = 2 -- 1D periodic; index = 3 -- 2D non-periodic; index = 4 -- 2D periodic; 

if index == 3 || index == 4 
    x = {[1,1], [1,2],[2,1],[2,2],[3,1],[1,3],[2,3],[3,2],[4,1],[1,4],[3,3],[3,4],[4,3]};
    x_labels = cellfun(@(c) sprintf('(%d, %d)', c(1), c(2)), x, 'UniformOutput', false);
    len = length(x) ;
else
    len = 6 ;
end

fable_nCNOT_metric = zeros(1,len) ;
fable_nG_metric = zeros(1,len) ;
fable_time = zeros(1,len) ;
bitble_nCNOT_metric = zeros(1,len) ;
bitble_nG_metric = zeros(1,len) ;
bitble_time = zeros(1,len) ;
bitble2_nCNOT_metric = zeros(1,len) ;
bitble2_nG_metric = zeros(1,len) ;
bitble2_time = zeros(1,len) ;
bitble3_nCNOT_metric = zeros(1,len) ;
bitble3_nG_metric = zeros(1,len) ;
bitble3_time = zeros(1,len) ;
if index == 1
    Qiskit_nCNOT_metric = [15, 83, 379, 1579, 6635, 26795] ;
    Qiskit_nG_metric = [11, 62, 308, 1312, 5610, 22827] ;
    for i = 1 : len
        Qiskit_nCNOT_metric(i) = Qiskit_nCNOT_metric(i) .* Laplace_mat(i).normalization_factor ;
        Qiskit_nG_metric(i) = Qiskit_nG_metric(i) .* Laplace_mat(i).normalization_factor ;
    end
elseif index == 2
    Qiskit_nCNOT_metric = [15, 83, 379, 1609, 6635, 26923] ;
    Qiskit_nG_metric = [11, 62, 310, 1346, 5609, 22827] ;
    for i = 1 : len
        Qiskit_nCNOT_metric(i) = Qiskit_nCNOT_metric(i) .* Laplace_mat(i).normalization_factor ;
        Qiskit_nG_metric(i) = Qiskit_nG_metric(i) .* Laplace_mat(i).normalization_factor ;
    end
elseif index == 3
    Qiskit_nCNOT_metric = [4 ,60 ,57, 377, 315 ,322 ,1611 ,1611 ,1484, 1498 , 6571, 26923, 26923] ;
    Qiskit_nG_metric = [2, 37, 35, 302, 237, 238, 1352, 1351, 1213, 1239, 5542, 22795, 22825] ;
    for i = 1 : len
        Qiskit_nCNOT_metric(i) = Qiskit_nCNOT_metric(i) .* Laplace_mat(i).normalization_factor ;
        Qiskit_nG_metric(i) = Qiskit_nG_metric(i) .* Laplace_mat(i).normalization_factor ;
    end
elseif index == 4
    Qiskit_nCNOT_metric = [8, 83, 83, 379, 379, 379, 1611, 1611, 1611, 1609, 6635, 26923, 26923] ;
    Qiskit_nG_metric = [4, 64, 62, 312, 306, 304, 1352, 1351, 1353, 1341, 5604, 22827, 22827] ;
    for i = 1 : len
        Qiskit_nCNOT_metric(i) = Qiskit_nCNOT_metric(i) .* Laplace_mat(i).normalization_factor ;
        Qiskit_nG_metric(i) = Qiskit_nG_metric(i) .* Laplace_mat(i).normalization_factor ;
    end
end


for i = 1:len
    fable_nCNOT_metric(i) = record_Laplace(index, i).fable_nCNOT_metric ;
    fable_nG_metric(i) = record_Laplace(index, i).fable_nG_metric ;
    fable_time(i) = record_Laplace(index, i).fable_time ;
    bitble_nCNOT_metric(i) = record_Laplace(index, i).bitble_nCNOT_metric ;
    bitble_nG_metric(i) = record_Laplace(index, i).bitble_nG_metric ;
    bitble_time(i) = record_Laplace(index, i).bitble_time ;
    bitble2_nCNOT_metric(i) = record_Laplace(index, i).bitble2_nCNOT_metric ;
    bitble2_nG_metric(i) = record_Laplace(index, i).bitble2_nG_metric ;
    bitble2_time(i) = record_Laplace(index, i).bitble2_time ;
    bitble3_nCNOT_metric(i) = record_Laplace(index, i).bitble3_nCNOT_metric ;
    bitble3_nG_metric(i) = record_Laplace(index, i).bitble3_nG_metric ;
    bitble3_time(i) = record_Laplace(index, i).bitble3_time ;
end



LineWidth = 0.5 ;
color_gray = [0.7 0.7 0.7];
color_shallow_blue = [0.8, 0.8, 1.0];
color_gray2 = [0.9 0.9 0.9] ;

semilogy(1:len,bitble_nCNOT_metric,':ro','LineWidth',LineWidth,...
    'MarkerSize',MarkerSize1,'MarkerFaceColor','r');
% hold on;
% semilogy(1:len,bitble_nG_metric,'--rs','LineWidth',LineWidth,...
%     'MarkerSize',MarkerSize1,'MarkerFaceColor','r');
hold on;
semilogy(1:len,bitble2_nCNOT_metric,':bo','LineWidth',LineWidth,...
    'MarkerSize',MarkerSize2,'MarkerFaceColor','b');
% hold on;
% semilogy(1:len,bitble2_nG_metric,':bs','LineWidth',LineWidth,...
%     'MarkerSize',MarkerSize2,'MarkerFaceColor','b');
hold on;
semilogy(1:len,bitble3_nCNOT_metric,':o','LineWidth',LineWidth,...
    'color', color_shallow_blue, 'MarkerSize',MarkerSize2,'MarkerFaceColor', color_shallow_blue);
% hold on;
% semilogy(1:len,bitble3_nG_metric,'--s','LineWidth',LineWidth,...
%     'color', color_shallow_blue, 'MarkerSize',MarkerSize2,'MarkerFaceColor', color_shallow_blue);
hold on ;
semilogy(1:len,fable_nCNOT_metric,':o','LineWidth',LineWidth,...
    'Color',color_gray,'MarkerSize',MarkerSize1,'MarkerFaceColor',color_gray);
% hold on ;
% semilogy(1:len,Qiskit_nCNOT_metric,':o','LineWidth',LineWidth,...
%     'Color','k','MarkerSize',MarkerSize1,'MarkerFaceColor',color_gray2);
% hold on;
% semilogy(1:len,fable_nG_metric,':s','LineWidth',LineWidth,...
%     'Color',color_gray,'MarkerSize',MarkerSize2,'MarkerFaceColor',color_gray);

% hold on;
% yyaxis right ;
% y = zeros(len,4);
% labels = cell(len,4); 
% for i = 1:len
%     y(i,:) = [record_Laplace(i).bitble_time, record_Laplace(i).bitble2_time, record_Laplace(i).bitble3_time, record_Laplace(i).fable_time ];
%     labels(i,:) = {num2str(record_Laplace(i).bitble_time,"%.1f"), num2str(record_Laplace(i).bitble2_time,"%.1f"), num2str(record_Laplace(i).bitble3_time,"%.1f"), num2str(record_Laplace(i).fable_time,"%.1f")};
% end
% 
% 
% bar_fig = bar(y, "GroupWidth", 0.92) ;
% % set(gca, 'XTickLabel',str, 'XTick',1:numel(str)) ;
% bar_fig(1).FaceColor = 'r' ;
% bar_fig(2).FaceColor = 'b' ;
% bar_fig(3).FaceColor = color_shallow_blue ;
% bar_fig(4).FaceColor = color_gray;
% 
% % set(gca,'ytick',[],'yticklabel',[]);
% 
% 
% for i=1:length(bar_fig)
%     xtips = bar_fig(i).XEndPoints;
%     xtips = xtips(1:len) ;
%     ytips = bar_fig(i).YEndPoints;
%     ytips = ytips(1:len) ;
%     sublabels = labels(:,i);
%     % t = text(xtips,ytips,sublabels,'HorizontalLlignment','center',...
%     %     'VerticalLlignment','bottom') ;
%     if i == 1
%         t = text(xtips+0.1,ytips + 0.2,sublabels,'HorizontalLlignment','center',...
%         'VerticalLlignment','bottom','FontWeight','bold') ;
%     else
%         t = text(xtips+0.1,ytips + 0.2,sublabels,'HorizontalLlignment','center',...
%             'VerticalLlignment','bottom') ;
%     end
%     set(t,'Rotation',90);
% end
% 
% ylabel('time(s)','FontWeight','normal');
% ylim([0,2.*max(fable_time)]) ;
% axis on;

h2 = legend(["BITBLE$^1$ CNOT","BITBLE$^2$ CNOT","BITBLE$^3$ CNOT","FLBLE CNOT"],'Location','northwest','Interpreter','latex','Orientation','vertical');
set(h2,'Box','off');
ax = gca;
set(gca, 'XTick', 1:len);
% set(gca, 'XTickLabel', x_labels);

% ax.YLxis(1).Color = 'k';
% ax.YLxis(2).Color = 'k';

% x_labels = {2,3,4,5,6,7};
xlim([0.9,len+.1]);
% set(gca, 'XTick', 2:len);
if index == 3 || index == 4 
    set(gca, 'XTickLabel', x_labels);
    xlabel("$(n_x, n_y)$ qubit",'Interpreter','latex');
else
    xlabel("$n$ qubit",'Interpreter','latex');
end

ylabel('number of gates $\times$ normalization factor','Interpreter','latex');


h2 = legend(["BITBLE$^1$ CNOT","BITBLE$^2$ CNOT","BITBLE$^3$ CNOT","FABLE CNOT"],'Interpreter','latex','Location','southeast');
set(h2,'Box','off');
