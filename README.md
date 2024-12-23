# BITBLE-SIABLE
Binary Tree Block-Encoding (BITBLE) and Single Ancilla Block-Encoding (SIABLE): Two simulation-friendly block-encoding circuits. These two block-encoding protocols inherit a circuit compression algorithm form [FABLE](https://github.com/QuantumComputingLab/fable), which have a lower subnormalization factors and require less ancillary qubits.

These two algorithms are bulit on top of [QCLAB](https://github.com/QuantumComputingLab/qclab), and SIABLE is also bulit on top of [cosine-sine decomposition](https://www.mathworks.com/matlabcentral/fileexchange/50402-cosine-sine-decomposition).


### QCLAB - MATLAB Implementation ###

In order to run the MATLAB implementation of BITBLE or SIABLE:

1. Install [QCLAB](https://github.com/QuantumComputingLab/qclab)
2. Install [cosine-sine decomposition](https://www.mathworks.com/matlabcentral/fileexchange/50402-cosine-sine-decomposition)
3. Compile `csd()` by running `make_csd.m` in the file named `iwoodsawyer-csd-a23bac9`
4. Add `iwoodsawyer-csd-a23bac9` and `QCLAB` files to the `bitblt_siable-qclab` file or clone all these three files into your MATLAB path.

After installation, BITBLE can be run for a target matrix `A` as either:

 ```
logging = true ; % logging of this algorithm
offset = 0 ;     % Qubit offset of this quantum circuit
[circuit, subnormalization_factor, info] = bitble( A, 'cutoff', 1e-4, logging, offset ) ;
[circuit, subnormalization_factor, info] = bitble( A, 'percentage', 80, logging, offset ) ;
```

SIABLE can be run with a similiar command: 

 ```
logging = true ; % logging of this algorithm
offset = 0 ;     % Qubit offset of this quantum circuit
[circuit, subnormalization_factor, info] = siable( A, 'cutoff', 1e-4, logging, offset ) ;
[circuit, subnormalization_factor, info] = siable( A, 'percentage', 80, logging, offset ) ;
```
    
The first option (`'cutoff'`) ignores coefficients smaller than `1e-4` in absolute value, the second option
(`'percentage'`) applies an 80% compression and only retains the 20% largest coefficients. The `'percentage'` and `logging` options are only available in the MATLAB version of BITBLE and SIABLE.

## Reference
