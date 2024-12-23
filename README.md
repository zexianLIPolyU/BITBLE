# BITBLE-SIABLE
Binary Tree Block-Encoding (BITBLE) and Single Ancilla Block-Encoding (SIABLE): Two simulation-friendly block-encoding circuits.These two block-encoding protocols have inherited a circuit compression algorithm form [FABLE](https://github.com/QuantumComputingLab/fable), which have a lower subnormalization factors and require less ancillary qubits.

These two algorithms are bulit on top of [QCLAB](https://github.com/QuantumComputingLab/qclab), and SIABLE are also bulit on top of [cosine-sine decomposition](https://www.mathworks.com/matlabcentral/fileexchange/50402-cosine-sine-decomposition).


### QCLAB - MATLAB Implementation ###

In order to run the MATLAB implementation of BITBLE or SIABLE:

1. Install [QCLAB](https://github.com/QuantumComputingLab/qclab)
2. Install [cosine-sine decomposition](https://www.mathworks.com/matlabcentral/fileexchange/50402-cosine-sine-decomposition)
3. Run `make_csd.m` in the file of `iwoodsawyer-csd-a23bac9`
4. Add `bitblt_siable-qclab` and `QCLAB` directory to the file `bitblt_siable-qclab` or clone these two file in your MATLAB path.

After installation, FABLE can be run for a target matrix `A` as either:

 ```
logging = true ; % logging of this algorithm
offset = 0 ;     % Qubit offset of this quantum circuit
[circuit, normalized_factor, info] = bitble( A, 'cutoff', 1e-4, logging, offset ) ;
[circuit, normalized_factor, info] = bitble( A, 'percentage', 80, logging, offset ) ;
```  
    
The first option (`'cutoff'`) ignores coefficients smaller than `1e-4` in absolute value, the second option
(`'percentage'`) applies an 80% compression and only retains the 20% largest coefficients. The `'percentage'` and `logging` options are only available in the MATLAB version of BITBLE and SIABLE.

## Reference
