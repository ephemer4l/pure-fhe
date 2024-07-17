# pure-fhe
Python implementations of the TFHE and CKKS fully homomorphic encryption schemes

## CKKS

Implementation in Python with Full RNS based on the [implementation](https://dspace.mit.edu/bitstream/handle/1721.1/129204/1227275316-MIT.pdf) by Saroja Erabelli (Available also on [github](https://github.com/sarojaerabelli/py-fhe)). Requires SymPy.

### TODO
- [ ] Get rid of SymPy
- [ ] Correct typing
- [ ] [Full-RNS](https://eprint.iacr.org/2018/931.pdf)
- [ ] [Efficient Bootstrapping](https://eprint.iacr.org/2020/1203)
- [ ] [SHARP](https://dl.acm.org/doi/abs/10.1145/3579371.3589053) accelerator
- [ ] [Taiyi](https://arxiv.org/abs/2403.10188) accelerator
- [ ] [HEAP](https://bu-icsg.github.io/publications/2024/fhe_parallelized_bootstrapping_isca_2024.pdf) accelerator

## TFHE

Implementation in Python based on [tfhe-py](https://github.com/pmuens/tfhe-py) by Philipp Muens. Requires Numpy.

## TODO

- [ ] Variable parameters using [lattice-estimator](https://github.com/malb/lattice-estimator)
- [ ] Correct typing
- [ ] [Public key encryption](https://eprint.iacr.org/2023/603.pdf)
