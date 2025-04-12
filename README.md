# Binary Tree Block-Encoding (BITBLE) and Single Ancilla Block-Encoding (SIABLE)
Binary Tree Block-Encoding (BITBLE) and Single Ancilla Block-Encoding (SIABLE): Two block-encoding quantum circuits for encoding general matrices. 
These two algorithms are bulit on top of [QCLAB](https://github.com/QuantumComputingLab/qclab) repository in MATLAB, and [mindquantum](https://www.mindspore.cn/mindquantum/docs/en/r0.6/index.html)/[PyQPanda](https://github.com/OriginQ/pyQPanda-Toturial/blob/master/source/index.rst) in Python

# 1. Binary Tree Block-Encoding (BITBLE)

Paper: https://arxiv.org/abs/2504.05624

Introduction:

BITBLE can compile a block-encoding very fast.

| Function      | Normalization Factor | Decoupled Methods                  | Ancilla Qubit Number for Encoding $2^n\times 2^n$ Matrix $A$    |
| -----------   | -----------          | -----------                        |  -----------                                                    |
| `bitble.m`    | $\Vert A\Vert_F$     | recursive demultiplexor            |  n                                                              |
| `bitble2.m`   | $\Vert A\Vert_F$     | permutative demultiplexor          |  n                                                              |
| `bitble3.m`   | $\mu_p(A^T)=\sqrt{\max_{i,j}\Vert A_{\cdot,j}\Vert_{2p}^{2p}\cdot \Vert A_{i,\cdot}\Vert_{2(1-p)}^{2(1-p)}}$  | recursive demultiplexor  |  n+2 |



## 1.1. Binary tree data structure state preparation - MATLAB Implementation ##



In order to run the MATLAB implementation of [state_preparation](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/state_preparation):
1. Down [state_preparation](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/state_preparation) repository.
2. Unzip it and add `state_preparation`  files into your MATLAB path.
    ```
    cd("state_preparation")
    ```
3. State preparation can be run by
     ```
     run("test_state_preparation.m")
     ```


## 1.2. Binary Tree Block-Encoding (BITBLE) - MATLAB Implementation ##

In order to run the MATLAB implementation of [BITBLE](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/bitble-qclab):

1. Down [BITBLE](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/bitble-qclab) repository.
2. Unzip it and add `QCLAB` files and `bitble.m` (or `bitble2.m` or `bitble3.m`) into your MATLAB path.
    ```
    cd("bitble-qclab")
    ```
3. BITBLE can be run for a target matrix `A` as either:
     ```
    run("test_bitble.m")
     ```
    or

    - `compr_type = 'cutoff'` and `compr_val = 1e-8` ignores coefficients smaller than `1e-8` in absolute value; `compr_type = 'percentage’`；
    - `compr_val = 80` applies an `80%` compression and only retains the `20%` largest coefficients.
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

    
    Output (The **Mac** system will display the content without any dislocation.):
    ```
         ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓                                                                                                                                                                                                                                                                                                                                                                                                   ┏━━┓ 
    q0: ━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●━⨯━━━━━●━━━━━━━━━━━━━●━━━━━━━━●━━━━━━●━┫RY┣━
         ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃                                                         ┃                                                         ┃                                                                                                                 ┃                                                                                                                   ┃ ┃     ┃             ┃        ┃      ┃ ┗━━┛ 
              ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃ ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┃ ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┃                                                                                                                 ┃                                                                                                                   ┃ ┃     ┃             ┃        ┃ ┏━━┓ ┃ ┏━━┓ 
    q1: ━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━⊕━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━╋━⨯━━━╋━━━━━━●━━━━━━╋━━━━━━●━⊕━┫RY┣━⊕━┫RY┣━
              ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃   ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃                                                         ┃                                                         ┃                                                         ┃                                                         ┃ ┃ ┃   ┃      ┃      ┃      ┃   ┗━━┛   ┗━━┛ 
              ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃        ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃ ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┃ ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┃ ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┃ ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┏━━┓   ┃ ┃ ┃   ┃ ┏━━┓ ┃ ┏━━┓ ┃ ┏━━┓ ┃ ┏━━┓          
    q2: ━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━╋━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━⊕━╋━╋━⨯━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━⊕━┫RY┣━━━━━━━━━━
              ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃        ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃   ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃   ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃   ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃ ┗━━┛ ┃   ┃ ┃ ┃   ┗━━┛   ┗━━┛   ┗━━┛   ┗━━┛          
              ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃        ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃        ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃        ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃        ┃      ┃      ┃      ┃      ┃      ┃      ┃      ┃   ┃ ┃ ┃                                      
    q3: ━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━━━━╋━━━━━━╋━━━━━━╋━━━━━━●━━━⨯━╋━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
              ┃      ┃      ┃             ┃      ┃      ┃             ┃      ┃      ┃             ┃      ┃      ┃               ┃      ┃      ┃             ┃      ┃      ┃             ┃      ┃      ┃             ┃      ┃      ┃               ┃      ┃      ┃             ┃      ┃      ┃               ┃      ┃      ┃             ┃      ┃      ┃               ┃      ┃      ┃             ┃      ┃      ┃            ┃ ┃                                      
              ┃      ┃      ┃             ┃      ┃      ┃             ┃      ┃      ┃             ┃      ┃      ┃               ┃      ┃      ┃             ┃      ┃      ┃             ┃      ┃      ┃             ┃      ┃      ┃               ┃      ┃      ┃             ┃      ┃      ┃               ┃      ┃      ┃             ┃      ┃      ┃               ┃      ┃      ┃             ┃      ┃      ┃            ┃ ┃                                      
    q4: ━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━━╋━━━━━━●━━━━━━╋━━━━━━━━━━━━⨯━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
              ┃             ┃             ┃             ┃             ┃             ┃             ┃             ┃               ┃             ┃             ┃             ┃             ┃             ┃             ┃             ┃               ┃             ┃             ┃             ┃               ┃             ┃             ┃             ┃               ┃             ┃             ┃             ┃              ┃                                      
              ┃             ┃             ┃             ┃             ┃             ┃             ┃             ┃               ┃             ┃             ┃             ┃             ┃             ┃             ┃             ┃               ┃             ┃             ┃             ┃               ┃             ┃             ┃             ┃               ┃             ┃             ┃             ┃              ┃                                      
    q5: ━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━●━━━━━━━━━━━━━━⨯━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
             
    ```
    and the number of quantum gates in block-encoding can be found in data `info`.



## 1.3. Binary Tree Block-Encoding (BITBLE) - PYTHON Implementation ##

Link: [https://github.com/ChunlinYangHEU/BITBLE_python](https://github.com/ChunlinYangHEU/BITBLE_python) .

-----------

# 2. Single Ancilla Block-Encoding (SIABLE):
Coming


## 2.1. Single Ancilla Block-Encoding (SIABLE) ##
SIABLE have optimal normalization factor $\Vert A\Vert_2$ with single ancilla qubit.

Coming soon.








