# Binary Tree Block-Encoding (BITBLE) and Single Ancilla Block-Encoding (SIABLE)
Binary Tree Block-Encoding (BITBLE) and Single Ancilla Block-Encoding (SIABLE): Two block-encoding quantum circuits for encoding general matrices. 
These two algorithms are bulit on top of [QCLAB](https://github.com/QuantumComputingLab/qclab) repository in MATLAB, and [mindquantum](https://github.com/mindspore-ai/mindquantum)/[PyQPanda](https://github.com/OriginQ/pyQPanda-Toturial/blob/master/source/index.rst) in Python

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

## 1.4. Performance Comparsion ##

Setting of ``bitble`` and [``fable``](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/test_fable) programe:

    ```
    offset = 0 ;
    logging = true ;
    compr_type = 'cutoff' ;
    circuit_sim = false ;
    ```

Time of compiling $2^{n}\times 2^{n}$ random matrices (compute the single parameters and CNOT logical on a 32GB RAM wins computer):

|  n            | `bitble1` | `bitble2` |`bitble3`| [`fable`](https://github.com/zexianLIPolyU/BITBLE-SIABLE_matlab/tree/main/test_fable)   |`Qiskit`|
|---------------|-----------|-----------|---------|-----------|--------|
| 5             | 0.017     | 0.010     | 0.014   | **0.008** | 6.97   |
| 6             | 0.031     | 0.035     | 0.054   | **0.028** | 46.3   |
| 7             | **0.104** | 0.113     | 0.211   | 0.112     | 221.2  |
| 8             | 0.642     | **0.464** | 0.872   | 0.484     | 669.7  |
| 9             | **1.683** | 1.980     | 3.637   | 1.984     | 767.4  |
| 10            | **7.137** | 8.368     | 15.21   | 8.604     | 7262   |
| 11            | **29.35** | 35.20     | 65.04   | 36.49     | -      |
| 12            | **120.4** | 147.1     | 269.5   | 148.7     | -      |
| 13            | **498.6** | 623.6     | 1125    | 625.5     | -      |
| 14            | **2103**  | 3119      | 5342    | 2831      | -      |

where `fable` in this experiment only compile the rotation's parameters, and `Qiskit` encoding $2^{n+1}\times 2^{n+1}$ random unitary by unitary synthesis.

References:

-FABLE:
[https://github.com/QuantumComputingLab/fable](https://github.com/QuantumComputingLab/fable)

-Qiskit's unitary synthesis:
[https://docs.quantum.ibm.com/api/qiskit/synthesis](https://docs.quantum.ibm.com/api/qiskit/synthesis).

-----------

# 2. Single Ancilla Block-Encoding (SIABLE):
Coming


## 2.1. Single Ancilla Block-Encoding (SIABLE) ##
SIABLE have optimal normalization factor $\Vert A\Vert_2$ with single ancilla qubit.

Coming soon.








