function [circuit, normalization_factor, info] = bitble2( A, compr_type, compr_val, logging, offset, circuit_sim )
% BITBLE2 -- Binary Tree Block Encodings using permutative multiplexed rotation with Less CNOT Gates.
%
% INPUT
% -----
% A:            matrix to be block encoded
% compr_type:   type of compression algorithm:
%                 * 'percentage' : compr_val between 0-100%, x% largest
%                                  coefficients retained
%                 * 'cutoff'     : compr_vall determines cutoff value, larger
%                                  coefficients retained
% compr_val:    input parameter for compression algorithm
% logging:      true/false, if true info will log information about compression
% offset:       starting position in the circuit. The default value for 'offset' is 0.
%
% OUTPUT
% ------    
% circuit:                      QCLAB circuit that block encodes A 
% normalization_factor:         subnormalization factor of this block-encoding
% info:                         struct containing some info on compression algorithm and circuit
%
% Copyright LI Zexian 
%
% Reference: BITBLE
%            Quantum Circuits for General Multiqubit Gates. 
%            Fast Approximate BLock Encodings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subnormalization
    normalization_factor = norm(A,'fro') ;
    N = size(A, 1) ;
    n = log2( N ) ;
    
    assert( N == 2^n ) ;
    assert( N == size(A,2) ) ;

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
        if logging
            info.datatype = 'real' ; 
        end
        
        % Compute rotation angles "Varphi"
        % Each column of "Varphi" is computed from each column of "A_magnitude" by Rotation-Y binary trees.
        Varphi = zeros(pow2(n)-1,pow2(n)+1);
        for col = 1:pow2(n)
            Varphi(:,col) = AngleCompute( 'RY', A(:,col), 1 );      
        end
        norm_Acol = sqrt(sum(A.^2,1))';
        norm_Acol = norm_Acol./norm(norm_Acol);
        Varphi( :,end ) = AngleCompute( 'RY', norm_Acol, 1 );
        
        if logging, info.vec.original.Varphi = Varphi ; end

        % transformed angles
        Varphi_global = Varphi( :,end );
        for i = 2:n
            row = pow2(i-1):pow2(i)-1;
            Varphi_global(row) = UniformlyRotationAngle(Varphi_global(row));
        end
        Varphi = Varphi(:,1:pow2(n)) ;
        a = zeros((pow2(n)+1)*(pow2(n)-1),1) ;
        for i = 1:n
            row = pow2(i-1):pow2(i)-1 ;
            b = Varphi(row,:)' ;
            a( (row(1) - 1)*(i~=1)*pow2(n) + 1 : row(end)*pow2(n) ) = UniformlyRotationAngle( b(:) ) ;
        end
        a(end-pow2(n)+2:end) = Varphi_global ;
        if logging
            info.vec.transformed = a ;
        end
        % Theshold vector according to compression criterion
        if strcmp( compr_type, 'percentage' )
            [ ~, sortIdx ] = sort( abs(a),'ascend') ;
            cutoff = floor( (compr_val/100.0) * N^2 - 1 ) ;
            a( sortIdx(1:cutoff) ) = 0 ;
            if logging
                info.vec.zeroed = cutoff ;
            end
        elseif strcmp( compr_type, 'cutoff' )
            if logging
                info.vec.zeroed = sum( abs(a) <= compr_val ) ;
            end
            a( abs(a) <= compr_val ) = 0 ;
            
        end
        Varphi = a;
        HatTheta = false;
        if logging
            info.vec.compressed.Varphi = Varphi ;
        end

    else % complex case
        if logging
            info.datatype = 'complex' ; 
        end
        A_norm = abs(A);
        A_phase = angle(A);
        Varphi = zeros( pow2(n)-1, pow2(n)+1 );
        HatTheta = zeros( pow2(n), pow2(n) );
        for col = 1:pow2(n)
            Varphi(:,col) = AngleCompute('RY', A_norm(:,col));
            HatTheta(:,col) = AngleCompute('RZ', A_phase(:,col));
        end
        norm_Acol = sqrt(sum(A_norm.^2,1))' ;
        norm_Acol = norm_Acol ./ norm(norm_Acol) ;
        Varphi( :,end ) = AngleCompute('RY',norm_Acol) ;
        if logging
            info.vec.original.Varphi = Varphi ;
            info.vec.original.HatTheta = HatTheta ;
        end

        % transformed angles
        Varphi_global = Varphi( :,end ) ;
        HatTheta_global = HatTheta(1,:) ;
        for i = 2:n
            row = pow2(i-1):pow2(i)-1;
            Varphi_global(row) = UniformlyRotationAngle(Varphi_global(row)) ;
        end
        Varphi = Varphi( :, 1:pow2(n) ) ;
        HatTheta( 1:pow2(n), : ) = HatTheta( [2:pow2(n),1], : ) ;
        a_Varphi = zeros( (pow2(n)+1)*(pow2(n)-1), 1 ) ;
        a_HatTheta = zeros( pow2(2*n), 1 ) ;
        for i = 1:n
            row = pow2(i-1):pow2(i)-1 ;
            b_Varphi = Varphi(row,:)' ;
            a_Varphi( (row(1) - 1)*(i~=1)*pow2(n) + 1 : row(end)*pow2(n) ) = UniformlyRotationAngle( b_Varphi(:) ) ; 
            b_HatThera = HatTheta(row,:)' ;
            a_HatTheta( (row(1) - 1)*(i~=1)*pow2(n) + 1 : row(end)*pow2(n) ) = UniformlyRotationAngle( b_HatThera(:) ) ;
        end
        a_Varphi(end-pow2(n)+2:end) = Varphi_global ;
        a_HatTheta(end-pow2(n)+1:end) = UniformlyRotationAngle( HatTheta_global(:) ) ;

        a = [ a_Varphi(:); a_HatTheta(:) ];
        if logging
            info.vec.transformed = a ;
        end
        % Theshold vector according to compression criterion
        if strcmp( compr_type, 'percentage' )
            [ ~, sortIdx ] = sort( abs(a),'ascend') ;
            cutoff = floor( (compr_val/100.0) * N^2 ) ;
            a( sortIdx(1:cutoff) ) = 0 ;
            if logging
                info.vec.zeroed = cutoff ;
            end
        elseif strcmp( compr_type, 'cutoff' )
            if logging
                info.vec.zeroed = sum( abs(a) <= compr_val ) ;
            end
            a( abs(a) <= compr_val ) = 0 ;
        end
        Varphi = a( 1 : (pow2(n)-1)*(pow2(n)+1) );
        HatTheta = a( (pow2(n)-1)*(pow2(n)+1)+1 : end );
        if logging
            info.vec.compressed.Varphi = Varphi ;
            info.vec.compressed.HatTheta = HatTheta ;
        end
    end
    % circuit
    [ circuit, info_circ ] = GenerateSingleControlQCircuit2( Varphi, HatTheta, offset, logging, circuit_sim ); 
    if logging
        info.circ = info_circ.circ ; 
    end
end % end of bitble2


function [ circuit, info ] = GenerateSingleControlQCircuit2( Varphi, Theta, offset, logging, circuit_sim )
% GenerateSingleControlQCircuit Generate a circuit which is consists of
% single-qubit controled gate and single-qubit rotation gate
% input:  Varphi   -- $2^n-1$-by-$2^n+1$ matrix, computed by Rotation-Y binary trees with "A_norm"
%         Theta    -- $2^n$-by-$2^n$ matrix, computed by Rotation-Z binary trees with "A_phase"
%         offset   -- default 0, the starting position this block-encoding
%         logging  -- true/false, if true info will log information about compression
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
    
    n = log( size(Varphi,1) + 1 )/(2*log( 2 )) ;
    Varphi_global = Varphi(end-pow2(n)+2:end) ;
    % Varphi = Varphi( :, 1:pow2(n) ) ;
    
    % The number of qubits
    NumQubits = 2*n;
    % The number of CNOT and Rotation-Z Rotation-Y
    if logging, nRZ = 0; nRY= 0; nCNOT = 0; nSWAP = 0; end 

    % circuit
    if circuit_sim
        circuit_PREP = qclab.QCircuit(NumQubits);
    else
        circuit_PREP = false ;
    end
    ctrl_index = n:(2*n-1) ;
    if iscomplex
        % Add the corresponding phase to each column 
        [ circuit_PREP, info_subcircuit ] = UniformRotation2( circuit_PREP, 'RZ', Theta(end-pow2(n)+1:end)', ctrl_index, 0, logging, circuit_sim ); 
        if logging
            nRZ = nRZ + info_subcircuit.nG ;
            nCNOT = nCNOT + info_subcircuit.nCNOT ;
        end
    end
    [ circuit_PREP, info_subcircuit ] = UniformRotation2( circuit_PREP, 'RY', Varphi(1:pow2(n))', ctrl_index, 0, logging, circuit_sim ) ;
    if logging
        nRY = nRY + info_subcircuit.nG ;
        nCNOT = nCNOT + info_subcircuit.nCNOT ; 
    end
    for k = 1:n-1
        para_row_index = pow2(k):pow2(k+1)-1 ;
        para_seq = Varphi( (para_row_index(1) - 1)*pow2(n) + 1 : para_row_index(end)*pow2(n) ) ;
        if any( para_seq ~= 0)
            ctrl_index = [0:k-1, n:(2*n-1)] ;
            [ circuit_PREP, info_subcircuit ] = UniformRotation2( circuit_PREP, 'RY', para_seq(:)', ctrl_index, k, logging, circuit_sim ) ; 
            if logging
                nRY = nRY + info_subcircuit.nG ;
                nCNOT = nCNOT + info_subcircuit.nCNOT ; 
            end
        end
    end
    if iscomplex
        [ circuit_PREP, info_subcircuit ] = UniformRotation2( circuit_PREP, 'RZ', Theta(1:pow2(n))', n:(2*n-1), 0, logging, circuit_sim ) ;
        if logging
            nRZ = nRZ + info_subcircuit.nG ;
            nCNOT = nCNOT + info_subcircuit.nCNOT ;
        end
        for k = 1:n-1
            para_row_index = pow2(k):pow2(k+1)-1 ;
            para_seq = Theta( (para_row_index(1) - 1)*pow2(n) + 1 : para_row_index(end)*pow2(n) ) ;
            if any( para_seq ~= 0)
                ctrl_index = [0:k-1, n:(2*n-1)] ;
                [ circuit_PREP, info_subcircuit ] = UniformRotation2( circuit_PREP, 'RZ', para_seq(:)', ctrl_index, k, logging, circuit_sim ) ; 
                if logging
                    nRZ = nRZ + info_subcircuit.nG ;
                    nCNOT = nCNOT + info_subcircuit.nCNOT ; 
                end
            end
        end
    end
    
    % Prepare the global norm
    if circuit_sim
        circuit_GLOBAL = qclab.QCircuit(NumQubits) ; 
        circuit_GLOBAL.push_back(qclab.qgates.RotationY(0,Varphi_global(1))) ; 
    else
        circuit_GLOBAL = false ;
    end
    if logging, nRY = nRY + 1 ; end 
    for k = 1 : n-1 
        para_row_index = pow2(k):pow2(k+1)-1 ;
        [circuit_GLOBAL, info_subcircuit ] = UniformRotation2( circuit_GLOBAL, 'RY', Varphi_global(para_row_index)', setdiff(0:n-1,k), k, logging, circuit_sim ) ; 
        if logging
            nRY = nRY + info_subcircuit.nG ; 
            nCNOT = nCNOT + info_subcircuit.nCNOT ; 
        end
    end
    
    if circuit_sim
        circuit = qclab.QCircuit(NumQubits, offset) ;
        circuit.push_back(circuit_PREP) ;
    else
        circuit = false ;
    end
    % Perform SWAP gates 
    if circuit_sim
        for swapq = 0 : NumQubits/2-1
            circuit.push_back(qclab.qgates.SWAP(swapq,NumQubits/2+swapq)) ;
        end
    end
    if logging, nSWAP = NumQubits/2 ; end 
    if circuit_sim
        circuit.push_back(circuit_GLOBAL.ctranspose) ;
    end
    if logging
        info.circ.nCNOT = nCNOT ;
        info.circ.nRY = nRY ;
        info.circ.nRZ = nRZ ;
        info.circ.nSWAP = nSWAP ;
    else
        info = false;
    end
end % end of GenerateSingleControlQCircuit2




%% UniformRotation22()


function [ circuit, info, parity_check ] = UniformRotation2( circuit, ctrl_type, para_seq, ctrl_index, targ_index, logging, circuit_sim )
% Input:    circuit     --  generated by qclab.QCircuit; 
%           ctrl_type   --  'RY' for Rotation-Y/ 'RZ' for Rotation-Z 
%           para_seq    --  parameter generated by Walsh-Hadamard transform 
%           ctrl_index  --  a vector contain the control 
%           logging     --  true/false, if true info will log information about compression 
% Output:   circuit     --  QCLAB circuit that block encodes A    
%           info        --  struct containing some info on compression algorithm and circuit

    n_count = log(size(para_seq, 2)) / log(2) ;
    if strcmp( ctrl_type, 'RY' )
        G = @qclab.qgates.RotationY ;
    elseif strcmp( ctrl_type, 'RZ' )
        G = @qclab.qgates.RotationZ ;
    end
    
    nG = 0 ; nCNOT = 0 ;
    i = 1;
    parity_check = int32(0);
    while i <= pow2(n_count)
        if any(para_seq(i) ~= 0)
            % Add CNOTs based on parity_check
            [ circuit, num_CNOT ] = make_CNOT(circuit, parity_check, ctrl_index, targ_index, circuit_sim ) ; 
            nCNOT = nCNOT + num_CNOT ; 

            if circuit_sim
                circuit.push_back( G(targ_index, para_seq(i)) ) ;
            end
            nG = nG + 1 ;
            ctrl = ctrl_pos( i, n_count ) ;
            % update parity check
            parity_check = bitset( int32(0), ctrl_index(ctrl) + 1, int32(1) ) ;
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
end % end of UniformRotation2

function ctrl = ctrl_pos(i, n)
    ctrl = n - log2(bitxor(grayCode(i-1),grayCode(i))) ;
    if i == pow2(n)
        ctrl = 1;
    end
end % end of ctrl_pos

function x = grayCode(x)
    x = bitxor(x,bitshift(x,-1));
end % end of grayCode


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



%% Uniformly Controlled Rotation Angle Computation
%  Reference: Quantum Circuits for General Multiqubit Gates. 2004.

function [thetat] = UniformlyRotationAngle(theta)
% Compute the uniformly controlled rotation
% thetat = (M^n)^(-1)*theta
    thetat = grayPermutation( sfwht( theta ) ) ;
end % end of UniformlyRotationAngle

function [ b ] = grayPermutation( a )
  k = log2( size(a, 1) ) ; 
  b = zeros( 2^k, 1 );
  for i = 0 : 2^k - 1
    b( i + 1 ) = a( grayCode( i ) + 1 );
  end
end % end of grayPermutation

function [ a ] = sfwht( a )
% Scaled fast Walsh-Hadamard transform
  k = log2(size(a, 1) ) ;
  for h = 1:k
    for i = 1:2^h:2^k
      for j = i:i+2^(h-1)-1
        x = a( j );
        y = a( j + 2^(h-1) );
        a( j ) = ( x + y ) / 2 ;
        a( j + 2^(h-1) ) = ( x - y ) / 2 ;
      end
    end
  end
end % end of sfwht



% function convertedAngles = AngleCompute( ctrl_type, NormOrPhase, is_real_leaves )
% 
%     if nargin <= 2
%         is_real_leaves = false ;
%     end
%     % Compute the Rotation Y angles "convertedAngle" by the Rotation Y binary tree
%     if strcmp( ctrl_type, 'RY' )
%         N = size(NormOrPhase, 1) ;
%         n = log(N) / log(2) ;
%         Varphi_seq = zeros( N-1, 1 ) ;
%         if is_real_leaves
%             [SumSquareRootAmplitude, Varphi] = positive_transform( NormOrPhase ) ;
%         else
%             [SumSquareRootAmplitude, Varphi] = AngleSearchBinTree( NormOrPhase );
%         end
%         Varphi_seq(pow2(n-1):pow2(n)-1) = Varphi ;
%         for i = n-1:-1:1
%             [SumSquareRootAmplitude, Varphi] = AngleSearchBinTree( SumSquareRootAmplitude ) ; 
%             Varphi_seq(pow2(i-1):pow2(i)-1) = Varphi ;
%         end
%         convertedAngles = Varphi_seq'.*2 ; 
%         % convertedAngles = mod( Varphi_seq'.*2, 4*pi ) ;
% 
%     % Compute the Rotation Z angles "convertedAngle" by the Rotation Z binary tree
%     elseif strcmp( ctrl_type, 'RZ' )
% 
%         InverseM = GenerateInversePhaseMatrix(length(NormOrPhase)) ;
%         convertedAngles = InverseM * NormOrPhase ; 
%         % convertedAngles = mod(convertedAngles, 4*pi) ;
% 
%     end
% end


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
        convertedAngles = NormOrPhase ;
    end
end

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
end

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
end


% %% Generate the matrix to compute Rotation Z angles
% function [InverseM] = GenerateInversePhaseMatrix(N)
% % Generate the phase calculation matrix M_k such that [theta] = M_k[hat_theta]
% % M_k is the matrix generated by bianry tree.
% % input: 'N' a interger of N = 2^n and N=k
% % 
% % Rotation Z binary tree:
% %                           -\hat{\theta_1}/2
% %                   /                                 \
% %         -\hat{\theta_2}/2                      +\hat{\theta_2}/2 
% %               /       \                         /             \
% %  -\hat{\theta_3}/2  +\hat{\theta_3}/2  -\hat{\theta_4}/2  +\hat{\theta_4}/2 
% %         =:\theta_1    =:\theta_2           =:\theta_3       =:\theta_4
% %      
% %  [\theta_1;\theta_2;\theta_3;\theta_4] = M_4 [\hat{\theta_1};\hat{\theta_2};\hat{\theta_3};\hat{\theta_4}]
% %  where M_{2k} = [kron(M_{k},[1;1], kron(eye(k),1/2.*[1;-1])].
% %  and InverseM_{2k} = [kron(InverseM_{k},[1,1]./2); kron(eye(2^(k-1)),[-1,1])].
% % 
% % ------------------------------------------------------------------------
%     InverseM = [-1, -1; -1, 1];
%     if N == 2
%         return;
%     end
%     k = 2;
%     while k ~= N
%         InverseM = [kron(InverseM,[1,1]./2); kron(eye(k),[-1,1])] ;
%         k = k * 2;
%     end
% end





