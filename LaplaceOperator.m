% 1D discretized Laplace operator:
% Lxx = [[ 2  -1   0  ...  * ];
%        [-1   2  -1  ... ...];
%        [ 0  ... ... ...  0 ];
%        [... ... -1   2  -1 ];
%        [ *  ...  0  -1   2 ]]
% where the entries ∗ in the lower-left and upper-right corner are
% either both equal to 0 for non-periodic boundary conditions,
% or both equal to −1 for periodic boundary conditions.

% In 2D, the discretized Laplace operator becomes the Kronecker sum of discretizations along the x- and y-directions:
%       L = Lxx ⊕ Lyy = Lxx ⊗ I + I ⊗ Lyy,
% which corresponds to a five-point stencil.

% Input: number of qubits:   n, 
%        dimision:           1 or 2, 
%        boundary_condition: 'on' or 'off' or ['on', 'on'] or ['on', 'off']
%                            or ['off', 'on'] or ['off', 'off']
% Output: 1D or 2D discretized Laplace operator
function L = LaplaceOperator(n, dimision, boundary_condition)
    % Check if the dimension is 1 or 2
    if dimision == 1
        % For 1D case, the boundary condition should be a logical value
        if ~islogical(boundary_condition)
            error('Boundary condition must be a logical value (true or false).');
        end
        L = oneDLapaceOperator(n, boundary_condition);
    elseif dimision == 2
        % For 2D case, the boundary condition should be a 2-element logical array
        if ~islogical(boundary_condition) || numel(boundary_condition) ~= 2
            error('Boundary condition must be a 2-element logical array with true or false.');
        end
        Lxx = oneDLapaceOperator(n(1), boundary_condition(1));
        Lyy = oneDLapaceOperator(n(2), boundary_condition(2));
        L = kron(Lxx, eye(pow2(n(2)))) + kron(eye(pow2(n(1))), Lyy);
    else
        error('Dimision should be 1 or 2.')
    end
end


% Input: number of qubits:   n, 
%        boundary_condition: 'on' or 'off'
% Output: 1D discretized Laplace operator
function L = oneDLapaceOperator(n, boundary_condition)
    % Initialize the 1D Laplacian operator
    L = eye(pow2(n));
    for i = 1:2^n-1
        L(i, i+1) = -1;
        L(i+1, i) = -1;
    end
    
    % Set boundary elements based on the boundary condition
    if boundary_condition
        L(1, end) = -1;     % Periodic boundary condition
        L(end, 1) = -1;
    else
        L(1, end) = 0;      % Non-periodic boundary condition
        L(end, 1) = 0;
    end
end