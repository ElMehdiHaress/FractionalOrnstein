# Fractional Ornstein–Uhlenbeck (fOU) — MATLAB/Octave

[![CI](https://github.com/ElMehdiHaress/FractionalOrnstein/actions/workflows/ci.yml/badge.svg)](https://github.com/ElMehdiHaress/FractionalOrnstein/actions/workflows/ci.yml)


Tools to simulate a 1-D fractional Brownian motion (fBm), build an fOU process, and reproduce estimator diagnostics (almost sure convergence + asymptotic distribution). Code runs in MATLAB or GNU Octave.

---

## Contents
├─ Almostsureconvergence_test.m # smoke/diagnostic test (lightweight)

├─ Asymptoticdistribution_test.m # smoke/diagnostic test (lightweight)

├─ eta_quantities.m # main helper for estimator (was: testest)

├─ fbm1d.m # Davies–Harte fBm sampler (1D)

├─ Newton2.m # small Newton solver used by tests

└─ (other .m files as needed)

## Requirements

- **MATLAB** R2020a+ (recommended), or  
- **GNU Octave** ≥ 8.0 (CI uses 8.4)

No external toolboxes are required.

---

## Quick start

### MATLAB
```matlab
% In MATLAB, from the repo root
addpath(genpath(pwd));           % add all .m files to the path

% Run the light tests
run('Almostsureconvergence_test.m');
run('Asymptoticdistribution_test.m');

% Example: compute estimator ingredients at level s with stepsize h
[result] = eta_quantities(8, 0.01);   % returns [H2, H3]
```

## Main functions

### `eta_quantities(s, h) -> [H2, H3]`

Computes the summary quantities used by the estimator from a simulated fOU path.

* `s` — level (uses `n = 2^s` coarse steps)
* `h` — macro step size
  Internally simulates fBm via `fbm1d`, builds the fOU with parameters set at the top of the file (edit `H`, `Theta`, `sigma` there).

### `fbm1d(H, N, T) -> B`

Generates a single-path 1D fBm of length `N+1` on `[0, T]` with Hurst `H` using the Davies–Harte method.

* `H ∈ (0,1)`, `N` time steps, `T` horizon
  Returns a column vector `B` with `N+1` points and `B(1)=0`.

### `Newton2(...)`

Minimal Newton solver used by the tests. See the function header for arguments.

---

## Performance tips

The tests are configured **light** for CI. To run heavier experiments locally:

* In `eta_quantities.m`, you can increase `N` (fine sub-steps), `s` (so `n=2^s`), or Monte-Carlo loops (if you add them) for higher fidelity.
* Long runs scale roughly linearly with `N * n`. Start small (`s=7..9`, moderate `N`) and grow as needed.

---

## Reproducibility

`fbm1d.m` uses `randn` internally. For reproducible runs, set the RNG before calls:

```matlab
rng(42, 'twister');   % MATLAB
% or in Octave:
rand("seed", 42); randn("seed", 42);
```

---

## CI

The repository includes a GitHub Actions workflow that runs the **light** tests on Octave to ensure the code compiles and the basic pipeline works. Heavier loops are commented out to keep CI fast.

---

## License

MIT. See `LICENSE`.

---

## Citation

If you use this code in academic work, please cite the repository:

```
@misc{FractionalOrnstein,
  author = {El Mehdi Haress, Yaozhong Hu},
  title  = {Fractional Ornstein–Uhlenbeck (fOU) MATLAB/Octave code},
  year   = {2025},
  url    = {https://github.com/ElMehdiHaress/FractionalOrnstein}
}
```

# Authors
Yaozhong Hu and El Mehdi Haress
 
