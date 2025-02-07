
clc;clear;close all

m = 5 ;
A = randn(7,8)+randn(7,8).*1i ; 
compr_val = 1e-8 ; 


%assert( size(A,1) == size(A,2) ) ; 
N = max( size(A) ) ;
n = ceil(log(N)/log(2)) ; 
% 0-padding
if min( size(A) ) ~= pow2(n)
    if N ~= pow2(n)
        A( N+1:pow2(n), N+1:pow2(n) ) = eye(pow2(n)-N) ; 
    else
        A( pow2(n), pow2(n) ) = 0 ;
    end
end
[UA, normalized_factor] = BlockEncoding(A) ; 



[C,S,V1,W1,V2,W2] = csd( UA(1:pow2(n),1:pow2(n)), UA(1:pow2(n),pow2(n)+1:pow2(n+1)), UA(pow2(n)+1:pow2(n+1),1:pow2(n)), UA(pow2(n)+1:pow2(n+1),pow2(n)+1:pow2(n+1)) );
C = diag(C) ;
index1 = abs(C) > compr_val ;

C( ~index1 ) = 0 ;
V_1 = zeros(size(V1)) ; W_1 = zeros(size(W1)) ;
V_1( :, index1 ) = V1( :, index1 ) ;
W_1( index1, : ) = W1( index1, : ) ;
if any(~index1) && V1(end, end) ~= 0
    V1( end, end ) = 0 ;
end
if any(~index1) && W1(end, end) ~= 0
    W1( end, end ) = 0 ;
end

S = diag(S) ;
index2 = abs(S) > compr_val ;
S( ~index2 ) = 0 ;

V_2 = zeros(size(V2)) ; W_2 = zeros(size(W2)) ;
V_2( :, index2 ) = V2( :, index2 ) ;
W_2( index2, : ) = W2( index2, : ) ;
if any(~index2) && V_2(end, end) ~= 0
    V_2( end, end ) = 0 ;
end
if any(~index2) && W_2(end, end) ~= 0
    W_2( end, end ) = 0 ;
end
V1 = V_1 ;
W1 = W_1 ;
V2 = V_2 ;
W2 = W_2 ;

C = diag(C) ;
S = diag(S) ;


p = pow2(n); q = pow2(n);
norm( (kron([1 0; 0 0],V1)+kron([0 0; 0 1],V2)) * [C,-S;S,C] * (kron([1 0; 0 0],W1')+kron([0 0; 0 1],W2')) -  UA )
norm([V1',zeros(p,q);zeros(q,p),V2'] * UA * [W1,zeros(p,q);zeros(q,p),W2] - [C,-S;S,C])

%% BlcokEncoding() and CosineSineDecomposition()

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

function [ C, V1, V2, W1, W2 ] = real_transform_csd( C, V1, V2, W1, W2 ) 
% transform the last round cosine-sine decomposition order to make diag([V1, V2])
% == diag([W1, W2]) == eye(2)

    C = (sign(V1) * sign(W1)) .* C ; 
    theta = acos(C) ;
    U = diag([1, V2]) * [ C, -sin(theta); sin(theta), C ] * diag([1, W2]) ; 
    
end

% 
%{
p=3;
q=3;
A = randn(p+q);
Q = orth(A);
[C,S] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q));
[C,S,U1,V1] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q));
[C,S,U1,V1,U2,V2] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q));


p=30;
q=8;
A = randn(p+q);
Q = orth(A);
[C,S] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q));
[C,S,U1,V1] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q));
[C,S,U1,V1,U2,V2] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q))


p=2;
q=2;
A = randn(p+q);
Q = orth(A);
[C,S] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q));
[C,S,U1,V1] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q));
[C,S,U1,V1,U2,V2] = csd(Q(1:p,1:p),Q(1:p,p+1:p+q),Q(p+1:p+q,1:p),Q(p+1:p+q,p+1:p+q));

norm(Q - blkdiag(U1,U2)*[C -S; S C]*blkdiag(V1,V2)')
norm(S.^2 + C.^2)-1

n=2;
A = randn(2*n);
Q = single(orth(A));
[C,S] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));
[C,S,U1,V1] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));
[C,S,U1,V1,U2,V2] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));

norm(Q - blkdiag(U1,U2)*[C -S; S C]*blkdiag(V1,V2)')
norm(S.^2 + C.^2)-1

n=2;
A = complex(randn(2*n),randn(2*n));
Q = orth(A);
[C,S] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));
[C,S,U1,V1] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));
[C,S,U1,V1,U2,V2] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));

norm(Q - blkdiag(U1,U2)*[C -S; S C]*blkdiag(V1,V2)')
norm(S.^2 + C.^2)-1

n=2;
A = single(complex(randn(2*n),randn(2*n)));
Q = single(orth(A));
[C,S] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));
[C,S,U1,V1] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));
[C,S,U1,V1,U2,V2] = csd(Q(1:n,1:n),Q(1:n,n+1:2*n),Q(n+1:2*n,1:n),Q(n+1:2*n,n+1:2*n));

norm(Q - blkdiag(U1,U2)*[C -S; S C]*blkdiag(V1,V2)')
norm(S.^2 + C.^2)-1

A = randn(10,5);
A = A*A';
[Q,R] = qr(A);
[C,S,U1,V1,U2,V2] = csd(Q(1:5,1:5),Q(1:5,6:10),Q(6:10,1:5),Q(6:10,6:10));
norm(Q - blkdiag(U1,U2)*[C -S; S C]*blkdiag(V1,V2)')
norm(S.^2 + C.^2)-1

A = complex(randn(10,5),randn(10,5));
A = A*A';
[Q,R] = qr(A);
[C,S,U1,V1,U2,V2] = csd(Q(1:5,1:5),Q(1:5,6:10),Q(6:10,1:5),Q(6:10,6:10));
norm(Q - blkdiag(U1,U2)*[C -S; S C]*blkdiag(V1,V2)')
norm(S.^2 + C.^2)-1

A = randn(10,5);
A = A*A';
[Q,R] = qr(single(A));
[C,S,U1,V1,U2,V2] = csd(Q(1:5,1:5),Q(1:5,6:10),Q(6:10,1:5),Q(6:10,6:10));
norm(Q - blkdiag(U1,U2)*[C -S; S C]*blkdiag(V1,V2)')
norm(S.^2 + C.^2)-1

A = complex(randn(10,5),randn(10,5));
A = A*A';
[Q,R] = qr(single(A));
[C,S,U1,V1,U2,V2] = csd(Q(1:5,1:5),Q(1:5,6:10),Q(6:10,1:5),Q(6:10,6:10));
norm(Q - blkdiag(U1,U2)*[C -S; S C]*blkdiag(V1,V2)')
norm(S.^2 + C.^2)-1
%}