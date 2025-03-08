# Single Ancilla Block-Encoding (SIABLE) and Binary Tree Block-Encoding (BITBLE) 
Single Ancilla Block-Encoding (SIABLE) and Binary Tree Block-Encoding (BITBLE): Two block-encoding circuits for encoding general matrices. 

| function      | normalization factor | Decoupled multiplexed rotations    | Ancilla qubit number for encoding $2^n\times 2^n$ matrix $A$  |
| -----------   | -----------          | -----------                        |  -----------                                        |
| `siable.m`    | $\Vert A\Vert_2$     | permutative multiplexed rotations  |  1                                                  |
| `bitble.m`    | $\Vert A\Vert_F$     | recursive multiplexed rotations    |  n                                                  |
| `bitble2.m`   | $\Vert A\Vert_F$     | permutative multiplexed rotations  |  n                                                  |
| `bitble3.m`   | $\mu_p(A^T)=\sqrt{\max_{i,j}\Vert A_{\cdot,j}\Vert_{2p}^{2p}\cdot \Vert A_{i,\cdot}\Vert_{2(1-p)}^{2(1-p)}}$  | recursive multiplexed rotations    |  n+2 |

These two algorithms are bulit on top of [QCLAB](https://github.com/QuantumComputingLab/qclab) repository in MATLAB, and SIABLE is also bulit on top of [cosine-sine decomposition](https://www.mathworks.com/matlabcentral/fileexchange/50402-cosine-sine-decomposition).
Python implementation can be found in [https://github.com/ChunlinYangHEU/BITBLE_python](https://github.com/ChunlinYangHEU/BITBLE_python).


## 1. Single Ancilla Block-Encoding (SIABLE) - MATLAB Implementation ##

SIABLE have optimal normalization factor $\Vert A\Vert_2$ with single ancilla qubit.

In order to run the MATLAB implementation of [SIABLE](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/siable-qclab):

1. Download [SIABLE](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/siable-qclab) repository.
2. Unzip it and add `QCLAB` and `iwoodsawyer-csd-a23bac9` files and `siable.m` into your MATLAB path.
3. Compile `csd()` by running `make_csd.m` in the file named `iwoodsawyer-csd-a23bac9` (The **Windows** version of MATLAB can compile conveniently).
    ```
    cd("iwoodsawyer-csd-a23bac9");
    ```
    ```
    run("make_csd.m");
    ```
    ```
    run("test_csd.m");
    ```
    If the screen output:
    ```
    ans =
    
        1.0000
    
    
    ans =
    
        1.0000
    ```
    then, `csd()` has been compilation.
4. After compilation, SIABLE can be tested by the following commands: 
    
     ```
    cd('..');
    run("test_siable.m");
     ```
    or 
    
    Define a matrix `A`:
    
     ```
    clc;clear;close all
    addpath("iwoodsawyer-csd-a23bac9"); % loading csd() 
    addpath("QCLAB");  % loading QCLAB

    n = 3 ;
    A = randn(pow2(n),pow2(n)) + 1j .* randn(pow2(n),pow2(n)) ;
    ```
    The first option (`'cutoff'`) ignores coefficients smaller than `1e-8` in absolute value,
    the second option (`'percentage'`) applies an `80%` compression and only retains the `20%` largest coefficients.
    ```
    fprintf("\nSIABLE Block Encoding");
    fprintf("\n------------------------------------------------------------ \n");
    fprintf("parameter computing... \n") ;
    
    offset = 0 ;
    logging = true ;
    compr_type = 'cutoff' ;%'percentage'; 
    compr_val = 1e-8 ;
    circuit_sim = true ;
    [circuit, normalization_factor, info] = siable( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
    % offset = 0 ;
    % logging = true ;
    % compr_type = 'percentage'; 
    % compr_val = 1e-8 ;
    % circuit_sim = true ;
    % [circuit, normalization_factor, info] = siable( A, compr_type, compr_val, logging, offset, circuit_sim ) ;
    ```
    Show the simulation result:
    ```
    fprintf("1.0001 * 2-norm of A = %f \n",1.0001 *norm(A,2)) ;
    fprintf("normalization_factor = %f \n",normalization_factor) ;
    M = circuit.matrix;
    fprintf("norm(normalization_factor.*M(1:pow2(n),1:pow2(n))-A) = %e \n",norm(normalization_factor.*M(1:pow2(n),1:pow2(n))-A)) ;
    if logging, info.circ; end
    circuit.draw();
    ```




## 2. Binary Tree Block-Encoding (BITBLE) - MATLAB Implementation ##

In order to run the MATLAB implementation of [BITBLE](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/bitble-qclab):

1. Down [BITBLE](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/bitble-qclab) repository.
2. Unzip it and add `QCLAB` files and `bitble.m` (or `bitble2.m` or `bitble3.m`) into your MATLAB path.

BITBLE can be run for a target matrix `A` as either:
 ```
run("test_bitble.m")
 ```
or
 ```
clc;clear;close all;
addpath( 'QCLAB' ) ;
N = pow2(3) ;
A = randn(N, N) ;
logging = true ; % logging of this algorithm
offset = 0 ;     % Qubit offset of this quantum circuit
circuit_sim = true ; % true/false, if false info will only compute the circuit's parameters; if true info will also simulate the quantum circuit
[circuit, subnormalization_factor, info] = bitble( A, 'cutoff', 1e-4, logging, offset, circuit_sim ) ;
% [circuit, subnormalization_factor, info] = bitble( A, 'percentage', 80, logging, offset, circuit_sim ) ;
circuit.draw() ;
info
```
The first option (`'cutoff'`) ignores coefficients smaller than `1e-4` in absolute value, the second option
(`'percentage'`) applies an 80% compression and only retains the 20% largest coefficients. The `'percentage'` and `logging` options are only available in the MATLAB version of BITBLE and SIABLE.

## 3. Binary Tree Block-Encoding (BITBLE) - PYTHON Implementation ##

Link: [https://github.com/zexianLIPolyU/BITBLE_python](https://github.com/zexianLIPolyU/BITBLE_python).




## Reference
