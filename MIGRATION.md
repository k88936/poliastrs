# Poliastro Rust MVP Migration - Executive Summary

## Objective Completed ✓
Comprehensive inventory of the Python poliastro API surface completed with a prioritized MVP migration roadmap for Rust.

**Deliverable**: `POLIASTRO_PYTHON_API_INVENTORY.md` (14 KB, complete analysis)

---

## Key Findings

### API Surface Statistics
- **Source Code**: 35,390 lines of Python across 67 files
- **Test Coverage**: 16,750 lines of tests across 41 files
- **Test Count**: ~200+ individual test cases
- **Core Modules**: 12 major functional domains

### Complexity Breakdown
| Domain | LOC | Complexity | Test Count |
|--------|-----|------------|-----------|
| Orbit/States | 1,430 | Medium | 85 |
| Propagation (core) | 1,509 | High | 15+ |
| Propagation (wrappers) | 1,500 | Low | - |
| Elements/Angles | 1,036 | Medium | 30+ |
| Maneuvers | 8,033 | Medium | 15+ |
| Plotting | 1,500+ | Medium | 40+ |
| Constants/Bodies | 1,690 | Low | 8 |
| **MVP Total (Tier 1)** | **~5,000** | **Low-Medium** | **~150** |

---

## Technical Considerations

### **Language Bindings**
- Start with pure Rust library
- Add PyO3 bindings for Python compatibility (later)
- Consider WASM for browser usage

### **Dependency Recommendations**
| Feature | Crate | Why |
|---------|-------|-----|
| Linear Algebra | `nalgebra` | Matrix ops for conversions |
| Root Finding | `brent-find` | Anomaly conversions (Newton) |
| Special Functions | `special` | Elliptic integrals for anomalies |
| Units | `uom` | Quantity tracking (orbital mechanics) |
| Time | `chrono` + custom | Julian date handling |

### **Known Challenges**
1. **Astropy Integration**: Python uses astropy for units/time; Rust needs replacement
2. **Numba JIT**: Python uses @jit for speed; Rust is natively fast
3. **Array Broadcasting**: NumPy style; Rust needs explicit loops or ndarray
4. **Quaternions**: Less needed in basic two-body; defer to Phase 2

---

## Test Migration TODO Checklist (Mandatory Full Parity)

Source of truth: `py-counterpart/tests/**/test_*.py`  
Requirement: every file below must be mapped to Rust tests and passing.

- [x] `py-counterpart/tests/test_bodies.py`
- [ ] `py-counterpart/tests/test_czml.py` // ignore by now
- [x] `py-counterpart/tests/test_ephem.py`
- [x] `py-counterpart/tests/test_examples.py`
- [x] `py-counterpart/tests/test_frames.py`
- [x] `py-counterpart/tests/test_hyper.py`
- [x] `py-counterpart/tests/test_iod.py`
- [x] `py-counterpart/tests/test_maneuver.py`
- [x] `py-counterpart/tests/test_sensors.py`
- [x] `py-counterpart/tests/test_spheroid_location.py`
- [x] `py-counterpart/tests/test_stumpff.py`
- [x] `py-counterpart/tests/test_util.py`
- [x] `py-counterpart/tests/tests_core/test_core_propagation.py`
- [x] `py-counterpart/tests/tests_core/test_core_util.py`
- [x] `py-counterpart/tests/tests_core/tests_threebody/test_cr3bp_quantities_calculations.py`
- [x] `py-counterpart/tests/tests_earth/test_earth_util.py`
- [x] `py-counterpart/tests/tests_earth/test_earthsatellite.py`
- [ ] `py-counterpart/tests/tests_earth/tests_atmosphere/test_coesa62.py` // ignored by now
- [x] `py-counterpart/tests/tests_earth/tests_atmosphere/test_coesa76.py`
- [ ] `py-counterpart/tests/tests_earth/tests_atmosphere/test_jacchia77.py`// ignored by now
- [ ] `py-counterpart/tests/tests_plotting/test_gabbard.py`
- [ ] `py-counterpart/tests/tests_plotting/test_misc.py`
- [ ] `py-counterpart/tests/tests_plotting/test_orbit_plotter.py`
- [ ] `py-counterpart/tests/tests_plotting/test_porkchop.py`
- [ ] `py-counterpart/tests/tests_plotting/test_tisserand.py`
- [x] `py-counterpart/tests/tests_spacecraft/test_spacecraft.py`
- [x] `py-counterpart/tests/tests_threebody/test_cr3bp_char_quant.py`
- [x] `py-counterpart/tests/tests_threebody/test_flybys.py`
- [x] `py-counterpart/tests/tests_threebody/test_restricted.py`
- [x] `py-counterpart/tests/tests_threebody/test_soi.py`
- [x] `py-counterpart/tests/tests_twobody/test_angles.py`
- [x] `py-counterpart/tests/tests_twobody/test_elements.py`
- [x] `py-counterpart/tests/tests_twobody/test_events.py`
- [x] `py-counterpart/tests/tests_twobody/test_mean_elements.py`
- [x] `py-counterpart/tests/tests_twobody/test_orbit.py`
- [x] `py-counterpart/tests/tests_twobody/test_perturbations.py`
- [x] `py-counterpart/tests/tests_twobody/test_propagation.py`
- [x] `py-counterpart/tests/tests_twobody/test_sampling.py`
- [x] `py-counterpart/tests/tests_twobody/test_states.py`
- [x] `py-counterpart/tests/tests_twobody/test_thrust.py`
