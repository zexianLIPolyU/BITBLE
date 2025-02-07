%CSD Cosine Sine Decomposition
% [C,S,U1,V1,U2,V2]=CSD(X11,X12,X21,X22) computes the CS decomposition of
% an M-by-M partitioned orthogonal/unitary matrix X:
% 
%                                 [  I  0  0 |  0  0  0 ]
%                                 [  0  C  0 |  0 -S  0 ]
%     [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]T
% X = [-----------] = [---------] [---------------------] [---------] .
%     [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
%                                 [  0  S  0 |  0  C  0 ]
%                                 [  0  0  I |  0  0  0 ]
% 
% X11 is P-by-Q. The orthogonal/unitary matrices U1, U2, V1, and V2 are
% P-by-P, (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S
% are R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
% which R = MIN(P,M-P,Q,M-Q).
%
% The mex calls the SORCSD/DORCSD/CUNCSD/ZUNCSD named LAPACK function.