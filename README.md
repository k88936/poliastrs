# Poliastrs: Rust Port of Poliastro

This is a Rust library for Astrodynamics, ported from the Python [poliastro](https://github.com/poliastro/poliastro) library.

## Features

- Orbital elements conversions
- Propagation (Keplerian, Cowell)
- Maneuvers (Hohmann, Bielliptic, Lambert)
- Plotting (using `plotters`)

## Examples

We have migrated several examples from the Python documentation to Rust. You can run them using `cargo run --example <name>`.

### Quickstart

Based on the [Python Quickstart Guide](https://docs.poliastro.space/en/stable/quickstart.html).
Demonstrates:
- Creating orbits from vectors and classical elements
- Propagation
- Hohmann transfer
- Solving Lambert's problem

```bash
cargo run --example quickstart
```

### Analyzing Hohmann Transfers

Based on "Studying Hohmann transfers" example.
Plots the Delta-V cost of Hohmann transfers as a function of target radius ratio.

```bash
cargo run --example hohmann_transfer
```
Output: `hohmann_transfer.png`

### Perturbed Propagation (Cowell's Method)

Based on "Studying non-keplerian orbits: perturbations" section of Quickstart.
Demonstrates propagating an orbit with a constant thrust perturbation using Cowell's method.

```bash
cargo run --example cowell_prop
```

## Migration Status

See `MIGRATION.md` for details on the migration progress.
