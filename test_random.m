%% test random matrix

clc;clear; close all
addpath( "QCLAB" ) ;
addpath( "bitble-qclab" ) ;
addpath( "test_fable" ) ;

offset = 0 ;
logging = 1 ;
compr_type = 'cutoff' ;%'percentage'; 
compr_val = 1e-8 ;
circuit_sim = false ;

record_random = struct() ;

for n = 2:14
    fprintf("n = %d:\n",n) ;
    disp("Fro-norm     Spectra-norm   ratio") ;
    N = pow2(n) ;
    A = randn(N, N) ;
    fprintf(" %e       %e       %e \n",norm(A,'fro'),norm(A,2),norm(A,'fro')/norm(A,2)) ;
    record_random(n-1).F_norm = norm(A,'fro');
    record_random(n-1).Spectral_norm = norm(A,2);
    record_random(n-1).norm_rotia = norm(A,'fro')/norm(A,2) ;
    record_random(n-1).max_element = max(max(abs(A))) ;
     %% Simulate the BITBLE quantum circuit 
    
    fprintf("\n\nBITBLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time = toc;
    record_random(n-1).dimension = n ;
    record_random(n-1).bitble_nRY = info.circ.nRY ;
    if ~isreal(A)
    record_random(n-1).bitble_nRZ = info.circ.nRZ ;
    end
    record_random(n-1).bitble_nCNOT = info.circ.nCNOT ;
    record_random(n-1).bitble_subnormalized_factor = subnormalization_factor ;
    record_random(n-1).bitble_time = time ;
    record_random(n-1).bitble_nG_merit = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_random(n-1).bitble_nCNOT_merit = subnormalization_factor .* info.circ.nCNOT ;

    fprintf( "number of RY + RZ gates * subnormalization factor = %e \n", subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ) ;
    fprintf( "number of CNOT gates * subnormalization factor = %e \n", subnormalization_factor .* info.circ.nCNOT ) ;
    
    %% Simulate the BITBLE2 quantum circuit 
    
    fprintf("\n\nBITBLE2 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble2( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time = toc;
    record_random(n-1).bitble2_nRY = info.circ.nRY ;
    if ~isreal(A)
    record_random(n-1).bitble2_nRZ = info.circ.nRZ ;
    end
    record_random(n-1).bitble2_nCNOT = info.circ.nCNOT ;
    record_random(n-1).bitble2_subnormalized_factor = subnormalization_factor ;
    record_random(n-1).bitble2_time = time ;
    record_random(n-1).bitble2_nG_merit = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_random(n-1).bitble2_nCNOT_merit = subnormalization_factor .* info.circ.nCNOT ;
    
    fprintf( "number of RY + RZ gates * subnormalization factor = %e \n", subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ) ;
    fprintf( "number of CNOT gates * subnormalization factor = %e \n", subnormalization_factor .* info.circ.nCNOT ) ;

    %% Simulate the BITBLE3 quantum circuit 
    
    fprintf("\n\nBITBLE3 Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, subnormalization_factor, info] = bitble3( A, 0.5, compr_type, compr_val, logging, offset, circuit_sim ) ;
    time = toc;
    record_random(n-1).bitble3_nRY = info.circ.nRY ;
    if ~isreal(A)
    record_random(n-1).bitble3_nRZ = info.circ.nRZ ;
    end
    record_random(n-1).bitble3_nCNOT = info.circ.nCNOT ;
    record_random(n-1).bitble3_subnormalized_factor = subnormalization_factor ;
    record_random(n-1).bitble3_time = time ;
    record_random(n-1).bitble3_nG_merit = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    record_random(n-1).bitble3_nCNOT_merit = subnormalization_factor .* info.circ.nCNOT ;
    
    fprintf( "number of RY + RZ gates * subnormalization factor = %e \n", subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ) ;
    fprintf( "number of CNOT gates * subnormalization factor = %e \n", subnormalization_factor .* info.circ.nCNOT ) ;
    
    %% Simulate the FABLE quantum circuit 
    
    fprintf("\n\nFABLE Block Encoding \n");
    fprintf("------------------------------------------------------------ \n");
    fprintf("parameter computing... \n");
    tic;
    [~, OA, alpha, info] = fable( A, compr_type, compr_val, logging, circuit_sim) ;
    time = toc;
    record_random(n-1).fable_nRY = info.circ.nRY ;
    if ~isreal(A)
    record_random(n-1).fable_nRZ = info.circ.nRZ ;
    end
    subnormalization_factor = max(max(abs(A))) .* pow2(n) ;
    record_random(n-1).fable_nCNOT = info.circ.nCNOT ;
    record_random(n-1).fable_subnormalized_factor = subnormalization_factor ;
    record_random(n-1).fable_time = time ;  
    if ~real(A)
        record_random(n-1).fable_nG_merit = subnormalization_factor .* (info.circ.nRY + info.circ.nRZ) ;
    else
        record_random(n-1).fable_nG_merit = subnormalization_factor .* info.circ.nRY ;
    end
    record_random(n-1).fable_nCNOT_merit = subnormalization_factor .* info.circ.nCNOT ;

    if logging
        info.circ
    end 
    nRG = info.circ.nRY ;
    if ~isreal(A)
        nRG = nRG + info.circ.nRZ ;
    end
    fprintf( "number of RY + RZ gates * subnormalization factor = %e \n", subnormalization_factor .* nRG ) ;
    fprintf( "number of CNOT gates * subnormalization factor = %e \n", subnormalization_factor .* info.circ.nCNOT ) ;

end
save("record_random.mat","record_random") ;
%}

%% Plot the time of computing parameters

clc;clear;close all;
load("record_random.mat") ;

LineWidth = 1 ;
MarkerSize1 = 3.5 ;
MarkerSize2 = 5.5 ;
color_gray = [0.7 0.7 0.7] ;
color_shallow_blue = [0.8, 0.8, 1.0];

dimension = size(record_random,1) ;
fable_time = zeros(1, dimension) ;
bitble_time = zeros(1, dimension) ;
bitble2_time = zeros(1, dimension) ;
bitble3_time = zeros(1, dimension) ;

for i = 1:dimension
    fable_time(i) = record_random(i).fable_time ;
    bitble_time(i) = record_random(i).bitble_time ;
    bitble2_time(i) = record_random(i).bitble2_time ;
    bitble3_time(i) = record_random(i).bitble3_time ;
end


figure;
plot(11:dimension+1,bitble_time(10:end),'-.sr','LineWidth',LineWidth,...
'MarkerSize',MarkerSize2,'MarkerFaceColor','r');
hold on;
plot(11:dimension+1,bitble2_time(10:end),'-.sb','LineWidth',LineWidth,...
'MarkerSize',MarkerSize2,'MarkerFaceColor','b');
hold on;
plot(11:dimension+1,bitble3_time(10:end),'-.s','LineWidth',LineWidth,...
'color',color_shallow_blue,'MarkerSize',MarkerSize2,'MarkerFaceColor', color_shallow_blue);
hold on;
plot(11:dimension+1,fable_time(10:end),'-.s','Color',color_gray,'LineWidth',LineWidth,...
'MarkerSize',MarkerSize1,'MarkerFaceColor',color_gray);


h2 = legend(["BITBLE$^1$","BITBLE$^2$","BITBLE$^3$","FABLE"],'Interpreter','latex','Location','north','Orientation','horizontal');
set(h2,'Box','off');
ax = gca;
set(gca, 'XTick', 11:dimension+1);
xlabel("$n$ qubits","Interpreter","latex");
ylabel('time(s)');

