function [circuit, subnormalization_factor, info] = siable( A, compr_type, compr_val, logging, offset, circuit_sim )
% SIABLE -- Single Ancilla Block Encodings.
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
% offset:       starting position in the circuit. The default value for 'offset' is 0
% circuit_sim:  true/false, if true info will simulate the circuit
%
% OUTPUT
% ------
% circuit:                      QCLAB circuit that block encodes A 
% subnormalization_factor:      subnormalization factor of this block-encoding
% info:                         struct containing some info on compression algorithm and circuit
%
% Copyright LI Zexian 2025.
%
% Reference: Optimal block encoding of general matrix
%            Quantum Circuits for General Multiqubit Gates
%            Fast Approximate BLock Encodings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if nargin <= 4
        offset = 0 ;
    end

    assert( size(A,1) == size(A,2) );
    N = size(A,1);
    n = ceil(log(N)/log(2));
    % 0-padding
    if N ~= pow2(n)
        A(N+1:pow2(n), N+1:pow2(n)) = eye(pow2(n)-N);
        % or A(pow2(n),pow2(n)) = 0;
    end
    % generate the
    [U, subnormalization_factor] = BlockEncoding(A);
    
    % Compute the rotation angles by CS decomposition
    Theta = zeros(pow2(n+1)-1, pow2(n)) ;
    C_record = zeros(pow2(n+1)-1, pow2(n)) ;
    Theta_index = zeros(1, pow2(n+1)-1) ;
    
    if isreal(A)
        if logging, info.datatype = 'real' ; end 
        for i = 1:n+1 
            for j = 1:pow2(i-1)
                zeta_index = zeta_i_j(i,j,n+1) ;
                for k = 1:pow2(i-1)
                    index_k = (k-1)*pow2(n+2-i)+1:k*pow2(n+2-i) ;
                    index_j = (j-1)*pow2(n+2-i)+1:j*pow2(n+2-i) ;
                    [ C,~,V1,W1,V2,W2 ] = csd( U( index_j(1:length(index_j)/2), index_k(1:length(index_k)/2) ),...
                            U( index_j(1:length(index_j)/2), index_k(length(index_k)/2+1:end) ),...
                            U( index_j(length(index_j)/2+1:end),  index_k(1:length(index_k)/2) ),...
                            U( index_j(length(index_j)/2+1:end),  index_k(length(index_k)/2+1:end) ) ) ;
                    C = diag(C) ; 
                    index_Theta = (k-1)*pow2(n+1-i)+1:k*pow2(n+1-i) ; 
                    C_record( zeta_index, index_Theta ) = C;
                    if i == n+1 
                        [ theta, V1, V2, W1, W2 ] = real_transform_csd( C, V1, V2, W1, W2 ) ; 
                        Theta( zeta_index, index_Theta ) = theta ;
                    else
                        Theta( zeta_index, index_Theta ) = acos(C') ;
                    end
                    U( index_j, index_k ) = [ V1, V2; W1', W2' ] ;
                    % Theta( zeta_index, index_Theta ) = acos(C') ;
                    Theta_index( zeta_index ) = i ;
                end
            end
        end
    else
        if logging, info.datatype = 'complex' ; end 
        U = complex(U);
        for i = 1:n+1 
            for j = 1:pow2(i-1)
                zeta_index = zeta_i_j(i,j,n+1);
                for k = 1:pow2(i-1)
                    index_k = (k-1)*pow2(n+2-i)+1:k*pow2(n+2-i);
                    index_j = (j-1)*pow2(n+2-i)+1:j*pow2(n+2-i);
                    [ C,~,V1,W1,V2,W2 ] = csd( complex(U( index_j(1:length(index_j)/2), index_k(1:length(index_k)/2) )),...
                            complex(U( index_j(1:length(index_j)/2), index_k(length(index_k)/2+1:end) )),...
                            complex(U( index_j(length(index_j)/2+1:end),  index_k(1:length(index_k)/2) )),...
                            complex(U( index_j(length(index_j)/2+1:end),  index_k(length(index_k)/2+1:end) )) );
                    C = diag(C);
                    U( index_j, index_k ) = [ V1, V2; W1', W2' ];
                    index_Theta = (k-1)*pow2(n+1-i)+1:k*pow2(n+1-i);
                    Theta( zeta_index, index_Theta ) = acos(C');
                    Theta_index( zeta_index ) = i;
                end
            end
        end
    end
    if logging
        info.vec.original.singular_value = C_record ;
        info.vec.original.Varphi = Theta ;
        info.vec.original.HatTheta = U ;
    end
    
    Theta = Theta.*2 ;
    Theta = UniformlyRotationAngle(Theta')';
    if compr_val ~= false
        Theta( abs(Theta) <= compr_val ) = 0 ;
    end
    InverseM = GenerateInversePhaseMatrix( size(U,1) ) ;
    U = InverseM * angle(U)' ; 
    if compr_val ~= false
        U( abs(U) <= compr_val) = 0 ; 
    end
    for i = 2:n+1
        U( pow2(i-1)+1:pow2(i), : ) = UniformlyRotationAngle( U( pow2(i-1)+1:pow2(i), : ) ) ;
    end
    a = [ Theta(:); U(:) ] ;

    if logging, info.vec.transformed = a ; end

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

    % transformed angles
    Theta = reshape( a( 1 : (pow2(n+1)-1)*pow2(n) ), pow2(n+1)-1, [] ) ;
    U = reshape( a( 1+(pow2(n+1)-1)*pow2(n) : end ), pow2(n+1), [] ) ; 
    if logging
        info.vec.compressed.Varphi = Theta ;
        info.vec.compressed.HatTheta = U ;
    end
    % generate the qcircuit
    % The number of qubits
    NumQubits = n + 1;

    if logging, nRY = 0 ; nRZ = 0 ; nCNOT = 0 ; end 

    % circuit
    if circuit_sim
        fprintf("Buliding circuit...");
        circuit = qclab.QCircuit(NumQubits, offset);
    else
        circuit = false ;
    end
    global_phase = 0;
    [ circuit, phase, info_subcircuit ] = ucrz(circuit, U(:,end), logging, circuit_sim); 
    if logging
        nRZ = nRZ + info_subcircuit.nRZ ;
        nCNOT = nCNOT + + info_subcircuit.nCNOT ;
    end
    global_phase = global_phase + phase ;

    %show the process 
    backNum = 0;
    fprintf(repmat('\b',1,backNum));         
    backNum = fprintf('%d/%d',1,pow2(NumQubits));
    
    for i = pow2(NumQubits)-1:-1:1
        [ circuit, info_subcircuit ] = double_ucry( circuit, Theta(i,:), Theta_index(i), logging, circuit_sim ) ; 
        if logging
            nRY = nRY + info_subcircuit.nG ;
            nCNOT = nCNOT + + info_subcircuit.nCNOT ; 
        end
        [ circuit, phase, info_subcircuit ] = ucrz( circuit, U(:,i), logging, circuit_sim ) ;
        if logging
            nRZ = nRZ + info_subcircuit.nRZ ;
            nCNOT = nCNOT + + info_subcircuit.nCNOT ;
        end
        global_phase = global_phase + phase;

        % show the process 
        fprintf(repmat('\b',1,backNum)) ;         
        backNum = fprintf('%d/%d',pow2(NumQubits)-i+1,pow2(NumQubits)) ; 

    end
    if circuit_sim
        circuit.push_back(qclab.qgates.RotationZ(0,global_phase));
    end
    if logging, nRZ = nRZ + 1 ; end 

    if logging
        info.circ.nCNOT = nCNOT ;
        info.circ.nRY = nRY ;
        info.circ.nRZ = nRZ ;
    else
        info = false;
    end
end




function [ circuit, info ] = UniformRotation( circuit, ctrl_type, para_seq, ctrl_index, targ_index, logging, circuit_sim )
% Input:    circuit     --  generated by qclab.QCircuit; 
%           ctrl_type   --  'RY' for Rotation-Y/ 'RZ' for Rotation-Z 
%           para_seq    --  parameter generated by Walsh-Hadamard transform 
%           ctrl_index  --  a vector contain the control 
%           logging     --  true/false, if true info will log information about compression 
% Output:   circuit     --  QCLAB circuit that block encodes A    
%           info        --  struct containing some info on compression algorithm and circuit

    n = log(size(para_seq, 2)) / log(2) ;
    if strcmp( ctrl_type, 'RY' )
        G = @qclab.qgates.RotationY ;
    elseif strcmp( ctrl_type, 'RZ' )
        G = @qclab.qgates.RotationZ ;
    end
    
    nG = 0 ; nCNOT = 0 ;
    i = 1;
    parity_check = int32(0);
    while i <= pow2(n)
        if any(para_seq(i) ~= 0)
            % Add CNOTs based on parity_check
            [ circuit, num_CNOT ] = make_CNOT(circuit, parity_check, ctrl_index, targ_index, circuit_sim ) ; 
            nCNOT = nCNOT + num_CNOT ; 
            
            if circuit_sim
                circuit.push_back( G(targ_index, para_seq(i)) ) ;
            end
            nG = nG + 1 ;
            ctrl = ctrl_pos( i, n ) ;
            % update parity check
            parity_check = bitset( int32(0), ctrl_index(ctrl) + 1, int32(1) ) ;
            i = i + 1;
        else
            % update parity check
            while i <= pow2(n) && all( para_seq(:,i) == 0 )
                ctrl = ctrl_pos(i, n);
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
end

function [ circuit, nCNOT ] = make_CNOT(circuit, parity_check, ctrl_index, targ_index, circuit_sim )
% Add CNOTs based on parity_check in the final
    if nargin == 3
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
end


function [circuit, global_phase, info ] = ucrz(circuit, theta, logging, circuit_sim)
%UCRZ   Uniformly Controlled Rotation along the Pauli-Z axis.
%   circuit = UCRZ(theta) constructs the quantum circuit that implements
%   the unitary
%
%                                 [exp(-i.*theta(1))                   ]
%   UCRz=exp(-.5i.*global_phase).*[                  .                 ]
%                                 [                   exp(-i.*theta(N))]
%
%   The length of the vector theta must be a power of 2
% -------------------------------------------------------------------------
% output: exp(-0.5i.*global_phase).*diag(circuit.matrix).' = exp(theta.*1i)

    if nargin == 2
        logging = false ;
    end
    N = numel(theta);
    n = ceil(log(N)/log(2));
    assert(N == pow2(n));

    if all(theta == 0)
        global_phase = 0 ;
        if logging
            info = struct() ;
            info.nRZ = 0 ;
            info.nCNOT = 0 ;
        else
            info = false ;
        end
        return;
    end
    
    % circuit
    if logging, nRZ = 0 ; nCNOT = 0; end
    if circuit_sim
        circuit.push_back(qclab.qgates.RotationZ(0,theta(2))); 
    end
    if logging, nRZ = nRZ + 1; end 
    for i = 2:n
        % theta(pow2(i-1)+1:pow2(i)) = grayPermutation(sfwht(theta(pow2(i-1)+1:pow2(i)))) ;
        [circuit, info]  = UniformRotation(circuit,'RZ',theta(pow2(i-1)+1:pow2(i))',0:i-2,i-1,logging, circuit_sim) ;
        if logging 
            nCNOT = nCNOT + info.nCNOT ;
            nRZ = nRZ + info.nG ;
        end
    end
    global_phase = theta(1);    
    if logging
        info = struct() ;
        info.nRZ = nRZ ;
        info.nCNOT = nCNOT ;
    else
        info = false ;
    end

end

function [circuit, info ] = double_ucry(circuit, theta, theta_index, logging, circuit_sim)

%UCRY   Uniformly Controlled Rotation along the Pauli-Y axis.
%   circuit = UCRY(theta) constructs the quantum circuit that implements
%   the unitary
%
%               [ c1        -s1          ]               [ c(1)       -s(1)                                       ]
%               [      .          .      ]               [     .           .                                      ]
%               [         cN         -sN ]               [      c(N/2)       -s(N/2)                              ]
%       UCRy1 = [ s1         c1          ]               [ s(1)        c(1)                                       ]
%               [      .          .      ]               [     .           .                                      ]
%               [         sN         cN  ]   ...         [      s(N/2)        c(N/2)                              ]
%                                                UCRy2 = [                          c(N/2+1)       -s(N/2+1)      ]
%                                                        [                                  .              .      ]
%                                                        [                                    c(N)           -s(N)]
%                                                        [                          s(N/2+1)       c(N/2+1)       ]
%                                                        [                                  .              .      ]
%                                                        [                                    s(N)            c(N)]...
% -------------------------------------------------------------------------
%               [ c-s         ]
%               [ s c         ]
%               [    c-s      ]
%  ... UCRyn =  [    s c      ]
%               [        .    ]
%               [          c-s]
%               [          s c]
% -------------------------------------------------------------------------

    if nargin == 3
        logging = false ;
    end
    N = length(theta);
    n = log2(N) + 1;
    
    if all( theta == 0 ) 
        if logging
            info = struct() ;
            info.nG = 0 ;
            info.nCNOT = 0 ;
        else
            info = false ;
        end
        return ;
    end
        
    % circuit
    ctrl_index = setdiff(0:n-1,theta_index-1) ;
    [circuit, info ] = UniformRotation( circuit,'RY',theta,ctrl_index,theta_index-1, logging, circuit_sim );

end


function ctrl = ctrl_pos(i, n)
    ctrl = n - log2(bitxor(grayCode(i-1),grayCode(i))) ;
    if i == pow2(n)
        ctrl = 1;
    end
end



%%

% BlcokEncoding()
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
end

% zeta_i_j(): return the order of Rotation Y 
function zeta_value = zeta_i_j(i,j,n)
    zeta_value = pow2(n-i)*(2*j-1);
end


function [ theta, V1, V2, W1, W2 ] = real_transform_csd( C, V_1, V_2, W_1, W_2 ) 
% transform the last round cosine-sine decomposition order to make diag([V1, V2])
% == diag([W1, W2]) == eye(2)

    C = sign(V_1) .* C ; 
    theta = acos(C) ;
    theta = angle( C + sign(V_2)*sin(theta)*1i ) ;
    V1 = sign(V_1) * V_1 ; V2 = sign(V_2) * V_2 ;
    W1 = W_1 ; 
    W2 = W_2 * sign(V_2)* sign(V_1) ; 
    
end

function [thetat] = UniformlyRotationAngle(theta)
% Compute the uniformly controlled rotation
% thetat = (M^n)^(-1)*theta

    % {
    N = size(theta,1);
    n = log(N)/log(2);
    % n must be a interger
    assert(mod(n,1) == 0);
    
    % Generate the Walsh-Hadamard matrix times permutation matrix 'M^n': M(i,j)^n = -1^BinaryTimesGray(i,j)
    % Bitwise multipulication 'BinaryTimesGray': BinaryTimesGray(i,j) = b(i-1) bitwise multipule g(j-1)
    % where 'b()' are the bianry code and 'g()' are the gray code.
    BinaryTimesGray = bitand((0:N-1)',grayCode(0:N-1));
    M = zeros(N,N);
    
    % Counting the parity of count '1' of bianry code of each element in BinaryTimesGray
    for i = 1:n
    M = bitxor(M,bitand(BinaryTimesGray,ones(N,N)));
    BinaryTimesGray = bitshift(BinaryTimesGray,-1);
    end
    M = -1 .* M + (M == 0);
    thetat = pow2(-n) .* (M' * theta);
    %}
    % thetat = grayPermutation( sfwht( theta ) ) ;
end

function x = grayCode(x)
    x = bitxor(x,bitshift(x,-1));
end

function [InverseM] = GenerateInversePhaseMatrix(N)
% Generate the phase calculation matrix M_k such that [theta] = M_k[hat_theta]
% M_k is the matrix generated by bianry tree.
% input: 'N' a interger of N = 2^n and N=k
% 
% Rotation Z binary tree:
%                           -\hat{\theta_1}/2
%                   /                                 \
%         -\hat{\theta_2}/2                      +\hat{\theta_2}/2 
%               /       \                         /             \
%  -\hat{\theta_3}/2  +\hat{\theta_3}/2  -\hat{\theta_4}/2  +\hat{\theta_4}/2 
%         =:\theta_1    =:\theta_2           =:\theta_3       =:\theta_4
%      
%  [\theta_0;\theta_1;\theta_2;\theta_3] = M_4 [\hat{\theta_0};\hat{\theta_1};\hat{\theta_2};\hat{\theta_3}]
%  where M_{2k} = [kron(M_{k},[1;1], kron(eye(k),1/2.*[1;-1])].
%  and InverseM_{2k} = [kron(InverseM_{k},[1,1]./2); kron(eye(2^(k-1)),[-1,1])].
% 
% ------------------------------------------------------------------------
    InverseM = [-1, -1; -1, 1];
    if N == 2
        return;
    end
    k = 2;
    while k ~= N
        % upper_InverseM = kron(InverseM,[1,1]./2);
        % lower_InverseM = kron(eye(k),[-1,1]);
        % InverseM = [upper_InverseM; lower_InverseM];
        InverseM = [kron(InverseM,[1,1]./2); kron(eye(k),[-1,1])] ;
        k = k * 2;
    end
end
