# many_body_certification

This header-based C++ library implements semidefinite-programming relaxations
for quantum spin systems. It computes lower bounds on ground-state energies and
upper and lower bounds on observables. The implementation supports translation,
sign, and lattice symmetries and is used by the model drivers in the sibling
`steady_state_phase_transitions` repository.

The code was used for [Mapping phase diagrams of quantum spin systems through
semidefinite-programming relaxations](https://arxiv.org/abs/2507.03137).

## Installation

See [INSTALL.md](INSTALL.md) for complete macOS and Linux setup instructions.
It covers the compiler and library requirements, MOSEK license setup, all
required environment variables, verification commands, and integration with
`steady_state_phase_transitions`.

The short version, after installing the dependencies, is:

```sh
export CXX=clang++
export MOSEK=/path/to/mosek/platform-directory
export EIGEN=/path/to/include/eigen3
export XTENSOR=/path/to/include
export XTL=/path/to/include

make -C cpp test_normal_form_fast
./cpp/bin/test_normal_form_fast
```

Do not edit machine-local dependency paths into the Makefile. Set them in the
shell or in a local shell startup file instead.

## Repository layout

- `cpp/include/`: reusable lattice, operator, Hamiltonian, RDM, and SDP headers
- `cpp/include/sos/`: sum-of-squares and momentum-sector SDP routines
- `cpp/src/`: model programs
- `cpp/tests/`: focused examples and tests
- `cpp/bin/`: locally built executables (generated artifacts)

## Building

Build from this repository's `cpp` directory through its Makefile. For example,
the portable normal-form test is:

```sh
make -C cpp test_normal_form_fast
./cpp/bin/test_normal_form_fast
```

Targets are built independently and executables are written to `cpp/bin/`.
There is no default `all` target. The test and example sources are the most
reliable reference for constructing new models; start with
`cpp/tests/simple_spin_chain.cpp` for a small spin-system example. Some older
targets use platform-specific linker flags; [INSTALL.md](INSTALL.md) gives the
appropriate MOSEK-linked smoke target for macOS and Linux.

## Cleaning

```sh
make -C cpp clean
```

This removes every file under `cpp/bin/`, so copy any executable you need to
keep before running it.
