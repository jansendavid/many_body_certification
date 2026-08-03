# Installation

## Supported setup

The project is C++20 source built with Make. It is routinely used on macOS and
Linux. Native Windows builds are not currently provided; on Windows, use a
Linux environment such as WSL.

Required software:

- GNU Make
- a C++20 compiler (`clang++` is the most consistently supported choice)
- Eigen 3
- xtensor and xtl
- MOSEK with its C++ Fusion API and a valid license

MOSEK is proprietary. Install it and configure an academic or commercial
license according to the MOSEK instructions before running solver-based
programs.

## 1. Keep the repositories together

For programs in `steady_state_phase_transitions` to find this library, use this
layout:

```text
work-directory/
├── many_body_certification/
└── steady_state_phase_transitions/
```

The two directories may be cloned separately; they only need to share the same
parent directory.

## 2. Install the open-source dependencies

On macOS with Homebrew:

```sh
brew install llvm eigen xtensor xtl
```

Apple's system `clang++` may also work. Homebrew LLVM is useful when a newer
C++20 compiler is required.

On Ubuntu or Debian:

```sh
sudo apt update
sudo apt install build-essential clang make libeigen3-dev libxtensor-dev libxtl-dev
```

Install MOSEK separately. The `MOSEK` value used below is its platform-specific
directory: it must directly contain `h/`, `bin/`, and `src/fusion_cxx/`.

## 3. Configure environment variables

The Makefiles deliberately take dependency locations from environment
variables. They have the following exact meanings:

| Variable | Must point to a directory containing |
| --- | --- |
| `MOSEK` | `h/fusion.h`, `bin/libmosek64.*`, and `bin/libfusion64.*` |
| `EIGEN` | `Eigen/` |
| `XTENSOR` | `xtensor/` |
| `XTL` | `xtl/` |
| `CXX` | C++ compiler executable, normally `clang++` |
| `SDPDIR` | parent directory containing `many_body_certification/` |

Example with explicit, machine-specific paths:

```sh
export CXX=clang++
export MOSEK=/path/to/mosek/platform-directory
export EIGEN=/path/to/include/eigen3
export XTENSOR=/path/to/include
export XTL=/path/to/include
export SDPDIR=/path/to/work-directory
```

For a standard Homebrew installation, the include variables can instead be
derived without hard-coding the Homebrew prefix:

```sh
export CXX="$(brew --prefix llvm)/bin/clang++"
export EIGEN="$(brew --prefix eigen)/include/eigen3"
export XTENSOR="$(brew --prefix xtensor)/include"
export XTL="$(brew --prefix xtl)/include"
```

Set `MOSEK` and `SDPDIR` separately because their locations depend on where
MOSEK and the repositories were installed. Put these exports in a local shell
startup file if they should persist. Do not commit local absolute paths.

## 4. Validate the paths

From the parent `work-directory`, run:

```sh
test -f "$MOSEK/h/fusion.h" || echo "MOSEK is incorrect"
test -d "$EIGEN/Eigen" || echo "EIGEN is incorrect"
test -d "$XTENSOR/xtensor" || echo "XTENSOR is incorrect"
test -d "$XTL/xtl" || echo "XTL is incorrect"
test -d "$SDPDIR/many_body_certification/cpp/include" || echo "SDPDIR is incorrect"
"$CXX" --version
```

No output from the five `test` commands means that their paths passed.

## 5. Build and run smoke tests

First verify code that does not invoke the solver:

```sh
make -C many_body_certification/cpp test_normal_form_fast
./many_body_certification/cpp/bin/test_normal_form_fast
```

The executable should print `fast normal-form tests passed`.

Then compile a MOSEK-linked target. On macOS, use:

```sh
make -C many_body_certification/cpp simple_spin_chain
```

On Linux, use:

```sh
make -C many_body_certification/cpp test_spins
```

A compiler error about `fusion.h` indicates an incorrect `MOSEK` path. A
linker error mentioning `mosek64` or `fusion64` indicates that `MOSEK/bin` does
not contain libraries for the current operating system and CPU. The ANNNI
example in step 6 performs an actual solve; a license error there means the
build succeeded but the MOSEK license still needs to be configured.

## 6. Verify the steady-state ANNNI integration

With both repositories under `SDPDIR`, build the regular and correlation entry
points:

```sh
make -C steady_state_phase_transitions/code annni_model_chain_TI
make -C steady_state_phase_transitions/code annni_model_chain_TI_corr
./steady_state_phase_transitions/code/bin/annni_model_chain_TI --help
```

Run the supplied one-point example, which records every parameter in
`run_parameters.json` and in the output CSV:

```sh
python3 steady_state_phase_transitions/code/run_annni_example.py \
  --output-dir /tmp/annni_example \
  --degree 1
```

See `steady_state_phase_transitions/code/ANNNI_USAGE.md` for the full CLI and
larger production runs.

## Platform notes

- Some older Makefile targets explicitly select either `g++` or `clang++` and
  have platform-specific linker flags. The ANNNI targets and primary examples
  are currently macOS-oriented. If an older target fails while the smoke tests
  pass, inspect that target's command in `cpp/Makefile` before changing the
  global dependency variables.
- Build only the named target you need; there is no default target that builds
  the whole project.
- `make -C many_body_certification/cpp clean` removes all files in `cpp/bin/`.
