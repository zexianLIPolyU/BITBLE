# Binary Tree Block-Encoding (BITBLE) and Single Ancilla Block-Encoding (SIABLE)
Binary Tree Block-Encoding (BITBLE) and Single Ancilla Block-Encoding (SIABLE): Two simulation-friendly block-encoding circuits. These two block-encoding protocols inherit a circuit compression algorithm form [FABLE](https://github.com/QuantumComputingLab/fable), which have a lower subnormalization factors and require less ancillary qubits.

These two algorithms are bulit on top of [QCLAB](https://github.com/QuantumComputingLab/qclab), and SIABLE is also bulit on top of [cosine-sine decomposition](https://www.mathworks.com/matlabcentral/fileexchange/50402-cosine-sine-decomposition).


## 1. Binary Tree Block-Encoding (BITBLE) - MATLAB Implementation ##

In order to run the MATLAB implementation of BITBLE:

Add `QCLAB` files and `bitble.m` (or `bitble2.m` or `bitble3.m`) into your MATLAB path.

### `bitble.m` : Binary tree block encoding with normalization factor of $||A||_F$ using recursive multiplexed rotations
### `bitble2.m`: Binary tree block encoding with normalization factor of $||A||_F$ using permutative multiplexed rotations
### `bitble3.m`: Binary tree block encoding with normalization factor of S_q(A) using recursive multiplexed rotations

BITBLE can be run for a target matrix `A` as either:


 ```
N = pow2(3);
A = randn(N, N);
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

## 2. Binary Tree Block-Encoding (BITBLE) - PYTHON Implementation ##

Link: [BITBLE_python](https://github.com/zexianLIPolyU/BITBLE_python).



## 3. Single Ancilla Block-Encoding (SIABLE) - MATLAB Implementation ##

In order to run the MATLAB implementation of SIABLE:

1. Add `QCLAB` and `iwoodsawyer-csd-a23bac9` files and `siable.m` into your MATLAB path.
2. Compile `csd()` by running `make_csd.m` in the file named `iwoodsawyer-csd-a23bac9`

After compilation, SIABLE can be run with a similiar command: 

 ```
N = pow2(3);
A = randn(N, N);
logging = true ; % logging of this algorithm
offset = 0 ;     % Qubit offset of this quantum circuit
circuit_sim = true ; % true/false, if false info will only compute the circuit's parameters; if true info will also simulate the quantum circuit
[circuit, normalization_factor, info] = siable( A, 'cutoff', 1e-4, logging, offset, circuit_sim ) ;
% [circuit, normalization_factor, info] = siable( A, 'percentage', 80, logging, offset, circuit_sim ) ;
circuit.draw();
info
```
The first option (`'cutoff'`) ignores coefficients smaller than `1e-4` in absolute value, the second option
(`'percentage'`) applies an 80% compression and only retains the 20% largest coefficients. The `'percentage'` and `logging` options are only available in the MATLAB version of BITBLE and SIABLE.

## Reference
