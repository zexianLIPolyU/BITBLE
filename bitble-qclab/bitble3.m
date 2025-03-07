function [ circuit, normalization_factor, info] = bitble3( A, p, compr_type, compr_val, logging, offset, circuit_sim )
% BITBLE3 -- Binary Tree Block Encodings using recursive multiplexed rotation with normalization factor as 
% $\sqrt{\mu_{2p}(A^T)\mu_{2(1-p)}(A)} = \sqrt{\max_k\Vert A(:,k)\Vert_{2*p}*\max_k\Vert A(k,:)\Vert_{2*(1-p)}}$
%
% INPUT
% -----
% A:            matrix to be block encoded
% p:            constant of normalization factor
% compr_type:   type of compression algorithm:
%                 * 'percentage' : compr_val between 0-100%, x% largest
%                                  coefficients retained
%                 * 'cutoff'     : compr_vall determines cutoff value, larger
%                                  coefficients retained
% compr_val:    input parameter for compression algorithm
% logging:      true/false, if true info will log information about compression
% offset:       starting position in the circuit. The default value for 'offset' is 0.
% circuit_sim:  true/false, if true info will simulate the circuit, elseif
%               false will only compute the single-qubit rotation angles
%
% OUTPUT
% ------
% circuit:                      QCLAB circuit that block encodes A 
% normalization_factor:         normalization factor of this block-encoding
% info:                         struct containing some info on compression algorithm and circuit
%
% Copyright LI Zexian, YANG Chunlin 2024.
%
% Reference: BITBLE
%            Quantum Circuits for General Multiqubit Gates. 
%            Fast Approximate BLock Encodings. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    assert( 0<=p && p<=1 ) ;
    dominator_chi_R = sqrt( S_qA( abs(A)', 2*p) ) ;
    dominator_chi_L = sqrt( S_qA( abs(A), 2*(1-p) ) ) ;
    % subnormalization
    normalization_factor = dominator_chi_R * dominator_chi_L ;
    N = size(A, 1) ;
    n = log2( N ) ;
    
    assert( N == 2^n ) ;
    assert( N == size(A, 2) ) ;
    
    if nargin <= 5
        circuit_sim = true ;
    end
    if nargin <= 4
        offset = 0 ;
    end
    if nargin <= 3
        logging = false ;
    end
    if nargin == 1
        compr_type = 'percentage';
        compr_val = 100 ;
    end

    % info struct
    if logging
        info = struct() ;
    else
        info = false ;
    end
    
    if isreal(A)
        if logging, info.datatype = 'real' ; end 

        A_norm = sign(A) .* abs(A).^p ./ sqrt( sum( abs(A).^(2*p) ) )  ; 
        A_norm2 = (abs(A)').^(1-p) ./ sqrt( sum( abs(A)'.^(2*(1-p)) ) ) ; 

        Varphi = zeros( pow2(n)-1, pow2(n) );
        Varphi2 = zeros( pow2(n)-1, pow2(n) );
        for col = 1:pow2(n)
            Varphi(:, col) = AngleCompute('RY', A_norm(:,col), true); 
            Varphi2(:, col) = AngleCompute('RY', A_norm2(:,col)); 
        end
        chi_R = acos( sqrt( sum( abs(A).^(2*p)) ) ./ dominator_chi_R ) * 2  ;
        chi_L = acos( sqrt( sum( abs(A').^(2*(1-p))) ) ./ dominator_chi_L ) * 2 ;
        
        if logging
            info.vec.original.Varphi.Varphi = Varphi ;
            info.vec.original.Varphi.Varphi2 = Varphi2 ;
            info.vec.original.chi.chi_R = chi_R ;
            info.vec.original.chi.chi_L = chi_L ;
        end

        % transformed angles
        for i = 2:n
            row = pow2(i-1):pow2(i)-1;
            Varphi(row,:) = UniformlyRotationAngle(Varphi(row,:)) ; 
            Varphi2(row,:) = UniformlyRotationAngle(Varphi2(row,:)) ;
        end
        Varphi = UniformlyRotationAngle( Varphi(:,1:pow2(n))' )';
        Varphi2 = UniformlyRotationAngle( Varphi2(:,1:pow2(n))' )';
        chi_R = UniformlyRotationAngle( chi_R' );
        chi_L = UniformlyRotationAngle( chi_L' );
        c = [ Varphi(:); Varphi2(:); chi_R(:) ; chi_L(:) ];
        if logging
            info.vec.transformed = c ;
        end

        % Theshold vector according to compression criterion
        if strcmp( compr_type, 'percentage' )
            [ ~, sortIdx ] = sort( abs(c),'ascend') ;
            cutoff = floor( (compr_val/100.0) * N^2 - 1 ) ;
            c( sortIdx(1:cutoff) ) = 0 ;
            if logging
                info.vec.zeroed = cutoff ;
            end
        elseif strcmp( compr_type, 'cutoff' )
            if logging
                info.vec.zeroed = sum( abs(c) <= compr_val ) ;
            end
            c( abs(c) <= compr_val ) = 0 ;
        end
        HatTheta = false;
        Varphi = reshape( c(1 : (pow2(n)-1)*pow2(n) ), pow2(n)-1, [] );
        Varphi2 = reshape( c((pow2(n)-1)*pow2(n)+1 : pow2(2*n+1)-pow2(n+1) ), pow2(n)-1, [] );
        chi_R = reshape( c(pow2(2*n+1)-pow2(n+1)+1 : pow2(2*n+1)-pow2(n) ), pow2(n), [] );
        chi_L = reshape( c(pow2(2*n+1)-pow2(n) + 1 : pow2(2*n+1) ), pow2(n), [] );
        if logging
            info.vec.compressed.Varphi.Varphi = Varphi ;
            info.vec.compressed.Varphi.Varphi2 = Varphi2 ;
            info.vec.compressed.chi.chi_R = chi_R ;
            info.vec.compressed.chi.chi_L = chi_L ;
        end

    else % complex case
        if logging, info.datatype = 'complex' ; end 
        A_phase = angle(A); 
        A_norm = abs(A).^p ./ sqrt( sum( abs(A).^(2*p) ) )  ; 
        A_norm2 = (abs(A)').^(1-p) ./ sqrt( sum( abs(A)'.^(2*(1-p)) ) ) ; 

        chi_R = acos( sqrt( sum( abs(A).^(2*p)) ) ./ dominator_chi_R ) * 2  ;
        chi_L = acos( sqrt( sum( abs(A').^(2*(1-p))) ) ./ dominator_chi_L ) * 2 ;
        Varphi = zeros( pow2(n)-1, pow2(n) ) ;
        Varphi2 = zeros( pow2(n)-1, pow2(n) ) ;
        HatTheta = zeros( pow2(n), pow2(n) ) ;
        for col = 1:pow2(n)
            Varphi(:, col) = AngleCompute('RY', A_norm(:,col)); 
            Varphi2(:, col) = AngleCompute('RY', A_norm2(:,col));
            HatTheta(:, col) = AngleCompute('RZ', A_phase(:,col));
        end
        
        if logging
            info.vec.original.Varphi.Varphi = Varphi ;
            info.vec.original.Varphi.Varphi2 = Varphi2 ;
            info.vec.original.chi.chi_R = chi_R ;
            info.vec.original.chi.chi_L = chi_L ;
            info.vec.original.HatTheta = HatTheta ; 
        end
        
        % transformed angles
        HatTheta(1:N,:) = HatTheta([2:N,1],:);
        for i = 2:n
            row = pow2(i-1):pow2(i)-1;
            Varphi(row,:) = UniformlyRotationAngle(Varphi(row,:)) ; 
            Varphi2(row,:) = UniformlyRotationAngle(Varphi2(row,:)) ;
            HatTheta(row,:) = UniformlyRotationAngle(HatTheta(row,:)) ;
        end
        
        Varphi = UniformlyRotationAngle( Varphi(:,1:pow2(n))' )' ;
        Varphi2 = UniformlyRotationAngle( Varphi2(:,1:pow2(n))' )' ;
        chi_R = UniformlyRotationAngle( chi_R' ) ;
        chi_L = UniformlyRotationAngle( chi_L' ) ;
        HatTheta = UniformlyRotationAngle( HatTheta' )';
        c = [ Varphi(:); Varphi2(:); chi_R(:) ; chi_L(:); HatTheta(:) ];
        
        if logging, info.vec.transformed = c ; end 

        % Theshold vector according to compression criterion
        if strcmp( compr_type, 'percentage' )
            [ ~, sortIdx ] = sort( abs(c),'ascend') ;
            cutoff = floor( (compr_val/100.0) * N^2 ) ;
            c( sortIdx(1:cutoff) ) = 0 ;
            if logging
                info.vec.zeroed = cutoff ;
            end
        elseif strcmp( compr_type, 'cutoff' )
            if logging
                info.vec.zeroed = sum( abs(c) <= compr_val ) ;
            end
            c( abs(c) <= compr_val ) = 0 ;
        end
        Varphi = reshape( c(1 : (pow2(n)-1)*pow2(n) ), pow2(n)-1, [] );
        Varphi2 = reshape( c((pow2(n)-1)*pow2(n)+1 : pow2(2*n+1)-pow2(n+1) ), pow2(n)-1, [] );
        chi_R = reshape( c(pow2(2*n+1)-pow2(n+1)+1 : pow2(2*n+1)-pow2(n) ), pow2(n), [] );
        chi_L = reshape( c(pow2(2*n+1)-pow2(n) + 1 : pow2(2*n+1) ), pow2(n), [] );
        HatTheta = reshape( c(pow2(2*n+1)+1 : pow2(2*n+1)+pow2(2*n) ), pow2(n), [] );
        if logging
            info.vec.compressed.Varphi.Varphi = Varphi ;
            info.vec.compressed.Varphi.Varphi2 = Varphi2 ;
            info.vec.compressed.chi.chi_R = chi_R ;
            info.vec.compressed.chi.chi_L = chi_L ;
        end
    end
    % circuit
    [ circuit, info_circ ] = GenerateSingleControlQCircuit3( Varphi, Varphi2, chi_R, chi_L, HatTheta, offset, logging, circuit_sim );
    if logging, info.circ = info_circ.circ; end
end % end of bitble3


function [ circuit, info ] = GenerateSingleControlQCircuit3(Varphi, Varphi2, chi_R, chi_L, Theta, offset, logging, circuit_sim)
% GenerateSingleControlQCircuit Generate a circuit which is consists of
% single-qubit controled gate and single-qubit rotation gate
% input:  Varphi   -- $2^n-1$-by-$2^n+1$ matrix, computed by Rotation-Y binary trees with "A_norm"
%         Varphi2  -- $2^n-1$-by-$2^n+1$ matrix, computed by Rotation-Y binary trees with "A_norm2"
%         chi_R(L) -- two $2^n$-by-$1$ matrix
%         Theta    -- $2^n$-by-$2^n$ matrix, computed by Rotation-Z binary trees with "A_phase"
%         offset   -- default 0, the starting position this block-encoding
%         logging  -- true/false, if true info will log information about compression
%         circuit_sim --  true/false, if true info will simulate the circuit, elseif false will only compute the single-qubit rotation angles
% output: circuit  -- QCLAB circuit that block encodes A    
%         info     -- struct containing some info on compression algorithm and circuit

    if nargin <= 2
        offset = 0; logging = false;
    elseif nargin == 3
        logging = false;
    end
    if ~Theta
        iscomplex = false;
    else
        iscomplex = true;
    end
    
    n = log( size(Varphi,1) + 1 )/log( 2 ) ;
    
    % The number of qubits
    NumQubits = 2*n + 2 ;
    % The number of CNOT and Rotation-Z Rotation-Y
    if logging, nRZ = 0; nRY= 0; nCNOT = 0; end 

    % circuit of U_R
    if circuit_sim
        circuit_U_R = qclab.QCircuit(NumQubits) ;
    else
        circuit_U_R = false ;
    end
    if iscomplex
        [circuit_U_R, info_subcircuit] = UniformRotation( circuit_U_R, 'RZ', Theta(pow2(n),:), int32((n+2):(2*n+1)), int32(0), logging, circuit_sim ) ; 
        if logging
            nRZ = nRZ + info_subcircuit.nG ;
            nCNOT = nCNOT + info_subcircuit.nCNOT ; 
        end
    end
    [circuit_U_R, info_subcircuit] = UniformRotation( circuit_U_R, 'RY', Varphi(1,:), int32((n+2):(2*n+1)), int32(0), logging, circuit_sim ) ; 
    if logging
        nRY = nRY + info_subcircuit.nG ;
        nCNOT = nCNOT + info_subcircuit.nCNOT ; 
    end
    %ctrl_index = 0:n-1;
    for k = 1:n-1
        para_row_index = pow2(k):pow2(k+1)-1;
        % Recursion of UniformRotation()
        [circuit_U_R, info_subcircuit ] = UniformRotation( circuit_U_R, 'RY', Varphi(para_row_index,:)', int32(0:(n-1)), int32(k), logging, circuit_sim ) ; 
        if logging
            nRY = nRY + info_subcircuit.nG ;
            nCNOT = nCNOT + info_subcircuit.nCNOT ; 
        end
    end
    if iscomplex
        [circuit_U_R, info_subcircuit] = UniformRotation( circuit_U_R, 'RZ', Theta(1,:), int32((n+2):(2*n+1)), int32(0), logging, circuit_sim ) ; 
        if logging
            nRZ = nRZ + info_subcircuit.nG ;
            nCNOT = nCNOT + info_subcircuit.nCNOT ; 
        end
        for k = 1:n-1
            para_row_index = pow2(k):pow2(k+1)-1;
            [circuit_U_R, info_subcircuit ] = UniformRotation( circuit_U_R, 'RZ', Theta(para_row_index,:)', int32(0:(n-1)), int32(k), logging, circuit_sim ) ; 
            if logging
                nRZ = nRZ + info_subcircuit.nG ;
                nCNOT = nCNOT + info_subcircuit.nCNOT ; 
            end
        end
    end
    [circuit_U_R, info_subcircuit ] = UniformRotation( circuit_U_R, 'RY', chi_R',  int32((n+2):(2*n+1)), int32(n), logging, circuit_sim ) ;
    if logging
        nRY = nRY + info_subcircuit.nG ;
        nCNOT = nCNOT + info_subcircuit.nCNOT ; 
    end

    % circuit of U_L
    if circuit_sim
        circuit_U_L = qclab.QCircuit(NumQubits) ;
    else
        circuit_U_L = false ;
    end
    [circuit_U_L, info_subcircuit] = UniformRotation( circuit_U_L, 'RY', Varphi2(1,:), int32((n+2):(2*n+1)), int32(0), logging, circuit_sim ) ; 
    if logging
        nRY = nRY + info_subcircuit.nG ;
        nCNOT = nCNOT + info_subcircuit.nCNOT ; 
    end
    for k = 1:n-1
        para_row_index = pow2(k):pow2(k+1)-1;
        % Recursion of UniformRotation()
        [circuit_U_L, info_subcircuit ] = UniformRotation( circuit_U_L, 'RY', Varphi2(para_row_index,:)', int32(0:(n-1)), int32(k), logging, circuit_sim ) ; 
        if logging
            nRY = nRY + info_subcircuit.nG ;
            nCNOT = nCNOT + info_subcircuit.nCNOT ; 
        end
    end
    [circuit_U_L, info_subcircuit ] = UniformRotation( circuit_U_L, 'RY', chi_L',  int32((n+2):(2*n+1)), int32(n+1), logging, circuit_sim ) ;
    if logging
        nRY = nRY + info_subcircuit.nG ;
        nCNOT = nCNOT + info_subcircuit.nCNOT ; 
    end
    
    if circuit_sim
        circuit = qclab.QCircuit(NumQubits, offset) ;
        circuit.push_back(circuit_U_R);
        % Perform SWAP gates 
        for swapq = 0 : n-1
            circuit.push_back(qclab.qgates.SWAP(swapq,n+2+swapq)) ;
        end
        circuit.push_back(circuit_U_L.ctranspose);
    else
        circuit = false ;
    end
    if logging
        info.circ.nCNOT = nCNOT ;
        info.circ.nRY = nRY ;
        info.circ.nRZ = nRZ ;
        info.circ.nSWAP = n ;
    else
        info = false;
    end
end % end of GenerateSingleControlQCircuit3


function [ circuit, info, parity_check ] = UniformRotation( circuit, ctrl_type, para_seq, ctrl_index, targ_index, logging, circuit_sim )
% Input:    circuit     --  generated by qclab.QCircuit; 
%           ctrl_type   --  'RY' for Rotation-Y/ 'RZ' for Rotation-Z 
%           para_seq    --  parameter generated by Walsh-Hadamard transform 
%           ctrl_index  --  a vector contain the index of control qubits 
%           targ_index  --  the index of controlled qubit
%           logging     --  true/false, if true info will log information about compression 
%           circuit_sim --  true/false, if true info will simulate the circuit, elseif false will only compute the single-qubit rotation angles
% Output:   circuit     --  QCLAB circuit that block encodes A    
%           info        --  struct containing some info on compression algorithm and circuit

    n_count = log(size(para_seq, 2)) / log(2) ;
    n = (circuit.nbQubits - 1) / 2 ;
    if strcmp( ctrl_type, 'RY' )
        G = @qclab.qgates.RotationY ;
    elseif strcmp( ctrl_type, 'RZ' )
        G = @qclab.qgates.RotationZ ;
    end
    
    nG = 0 ; nCNOT = 0 ;
    i = 1;
    parity_check = int32(0);
    while i <= pow2(n_count)
        if any(para_seq(:,i) ~= 0)
            % Add CNOTs based on parity_check
            [ circuit, num_CNOT ] = make_CNOT( circuit, parity_check, ctrl_index, targ_index, circuit_sim ) ;
            nCNOT = nCNOT + num_CNOT ;
            % Reset parity check
            parity_check = int32(0) ;
            if size(para_seq, 1) == 1
                if circuit_sim
                    circuit.push_back( G(targ_index, para_seq(i)) ) ; 
                end
                nG = nG + 1 ;
            else
                [ circuit, info_subcircuit, parity_check ] = UniformRotation( circuit, ctrl_type, para_seq(:,i)', n + 1 + ctrl_index, targ_index, logging, circuit_sim );
                if logging 
                    nG = nG + info_subcircuit.nG ;
                    nCNOT = nCNOT + info_subcircuit.nCNOT ; 
                end
            end
            ctrl = ctrl_pos( i, n_count ) ;
            % update parity check
            parity_check = bitset( parity_check, ctrl_index(ctrl) + 1, int32(1) ) ;
            i = i + 1;
        else
            % update parity check
            while i <= pow2(n_count) && all( para_seq(:,i) == 0 )
                ctrl = ctrl_pos(i, n_count);
                if bitget( parity_check, ctrl_index(ctrl) + 1 ) 
                    parity_check = bitset( parity_check, ctrl_index(ctrl) + 1, int32(0) ) ;
                else
                    parity_check = bitset( parity_check, ctrl_index(ctrl) + 1, int32(1) ) ;
                end
                i = i + 1;
            end
        end
    end
    % Add CNOTs based on parity_check in the final
    [ circuit, num_CNOT ] = make_CNOT(circuit, parity_check, ctrl_index, targ_index, circuit_sim ) ;
    nCNOT = nCNOT + num_CNOT ;
    
    if logging
        info = struct() ;
        info.nG = nG ;
        info.nCNOT = nCNOT ;
    else
        info = false ;
    end
end % end of UniformRotation

function ctrl = ctrl_pos(i, n)
    ctrl = n - log2(bitxor(grayCode(i-1),grayCode(i))) ;
    if i == pow2(n)
        ctrl = 1;
    end
end % end of ctrl_pos

function [ circuit, nCNOT ] = make_CNOT(circuit, parity_check, ctrl_index, targ_index, circuit_sim )
% Add CNOTs based on parity_check in the final
    if nargin <= 3
        targ_index = 0 ;
    end
    nCNOT = 0 ; 
    if parity_check ~= int32(0)
        for j = 1 : numel(ctrl_index)
            if bitget( parity_check, ctrl_index(j) + 1 )
                if circuit_sim
                    circuit.push_back( qclab.qgates.CNOT( ctrl_index(j), targ_index ) ) ; 
                end
                nCNOT = nCNOT + 1 ;
            end
        end
    end
end % end of make_CNOT

%% Anglecompute

function SqA = S_qA(A,q)
    % the qth power of the maximum q-norm of any row 
    N = size(A, 1) ;
    row_value = zeros(N, 1) ;
    for k = 1 : N
        row_value(k,:) = sum( abs(A(k,:)).^q ) ;
    end
    SqA = max(row_value) ;
end % end of S_qA


function convertedAngles = AngleCompute( ctrl_type, NormOrPhase, is_real_leaves )
    
    if nargin <= 2
        is_real_leaves = false ;
    end

    N = size(NormOrPhase, 1) ;
    n = log(N) / log(2) ;
    % Compute the Rotation Y angles "convertedAngle" by the Rotation Y binary tree
    if strcmp( ctrl_type, 'RY' )
        
        Varphi_seq = zeros( N-1, 1 ) ;
        if is_real_leaves
            [SumSquareRootAmplitude, Varphi] = positive_transform( NormOrPhase ) ;
        else
            [SumSquareRootAmplitude, Varphi] = AngleSearchBinTree( NormOrPhase );
        end
        Varphi_seq(pow2(n-1):pow2(n)-1) = Varphi ;
        for i = n-1:-1:1
            [SumSquareRootAmplitude, Varphi] = AngleSearchBinTree( SumSquareRootAmplitude ) ; 
            Varphi_seq(pow2(i-1):pow2(i)-1) = Varphi ;
        end
        convertedAngles = mod( Varphi_seq'.*2, 4*pi ) ;

    % Compute the Rotation Z angles "convertedAngle" by the Rotation Z binary tree
    elseif strcmp( ctrl_type, 'RZ' )
        for i = 1 : n
            len = pow2(n-i) ;
            a = NormOrPhase( 1 : 2*len ) ;
            for j = 1 : len  
                NormOrPhase(j) = (a(2*j-1) + a(2*j))./2 ;
                NormOrPhase(len+j) = - a(2*j-1) + a(2*j) ;
            end
        end
        NormOrPhase(1) = -a(1)-a(2) ;
        convertedAngles = NormOrPhase;
    end
end % end of AngleCompute

%% compute Rotation Y angles from the leaves in a binary tree and updating the the leaves
function [SumSquareRootAmplitude, Varphi_list] = AngleSearchBinTree( amplitude )
% Compute the rotation theta list from the leaves
% Input: a real vector "amplitude"  (Note that the length of "amplitude" have to be 2^n)
% 
% Rotation Y binary tree:
%                                        1
%                       /                                       \
%                *cos(Varphi_1)                             *sin(Varphi_1) 
%               /           \                              /           \
%     *cos(Varphi_2)        *sin(Varphi_2)        *cos(Varphi_3)         *sin(Varphi_3)
%     =:amplitude_norm(1)  =:amplitude_norm(2)  =:amplitude_norm(3)   =:amplitude_norm(4)
% ------------------------------------------------------------------------
    lengthAmplitude = length(amplitude);    
    if lengthAmplitude == 2
        % whose childern is leaf note
        if all(amplitude == 0)
            SumSquareRootAmplitude = 0;
            Varphi_list = 0;
        else
            SumSquareRootAmplitude = sqrt(amplitude(1)^2+amplitude(2)^2) ;
            Varphi_list = acos(amplitude(1)/SumSquareRootAmplitude) ;
        end
    else
        [SumSquareAmplitude1, Theta_list1] = AngleSearchBinTree( amplitude(1:lengthAmplitude/2) );
        [SumSquareAmplitude2, Theta_list2] = AngleSearchBinTree( amplitude(lengthAmplitude/2+1:lengthAmplitude) );
        SumSquareRootAmplitude = [SumSquareAmplitude1,SumSquareAmplitude2];
        Varphi_list = [Theta_list1, Theta_list2];
    end
end % end of AngleSearchBinTree

function [SumSquareAmplitude, Varphi] = positive_transform( SignedAmplitude )
% compute the leaves nodes (with signed number) on Rotation-Y tree in the real case
    N = size( SignedAmplitude, 1 );
    Varphi = zeros(1, N/2);
    SumSquareAmplitude = zeros(1, N/2);
    for i = 1:N/2
        if all( SignedAmplitude( 2*i-1:2*i, 1 ) == 0 )
            SumSquareAmplitude(i) = 0;
            Varphi(i) = 0;
        else
            num = SignedAmplitude( 2*i-1:2*i, 1 ) ;
            SumSquareAmplitude(i) = norm( num, 2 ) ;
            num = num ./ SumSquareAmplitude(i) ;
            complex_num = num(1) + num(2).*1j;
            Varphi(i) = mod( angle(complex_num), 2*pi ) ;
        end
    end
end % end of positive_transform


%% Uniformly Controlled Rotation Angle Computation
%  Reference: Quantum Circuits for General Multiqubit Gates. 2004.

function [thetat] = UniformlyRotationAngle(theta)
% Compute the uniformly controlled rotation
% thetat = (M^n)^(-1)*theta
    
    thetat = grayPermutation( sfwht( theta ) );
end % end of UniformlyRotationAngle

function [ b ] = grayPermutation( a )
  k = log2( size(a, 1) ) ; 
  b = zeros( size(a) );
  for i = 0 : 2^k - 1
    b( i + 1, : ) = a( grayCode( i ) + 1, : );
  end
end % end of grayPermutation

function [ a ] = sfwht( a )
% Scaled fast Walsh-Hadamard transform
  k = log2(size(a, 1) ) ;
  for h = 1:k
    for i = 1:2^h:2^k
      for j = i:i+2^(h-1)-1
        x = a( j ,: );
        y = a( j + 2^(h-1) ,: );
        a( j ,: ) = ( x + y ) / 2 ;
        a( j + 2^(h-1) ,: ) = ( x - y ) / 2 ;
      end
    end
  end
end % end of sfwht



function x = grayCode(x)
    x = bitxor(x,bitshift(x,-1));
end % end of grayCode





