# Poliastro Python API Inventory and MVP Rust Migration Shortlist

## Project Overview
- **Total Source Code**: 35,390 lines of Python
- **Total Tests**: 16,750 lines of Python
- **Total Files**: 67 Python files across 22 directories
- **Test Files**: 41 test files with comprehensive coverage

## Core Package Structure

### 1. **Bodies & Constants** (`bodies.py`, `constants/`)
- **Location**: `/py-counterpart/src/poliastro/bodies.py` (373 lines)
- **Location**: `/py-counterpart/src/poliastro/constants/` (1,317 lines total)
- **Purpose**: Defines celestial bodies (Sun, Earth, Moon, planets, moons) and physical constants
- **Key Classes**: `Body`, `SolarSystemPlanet`
- **Core Data**: Gravitational parameters, radii, rotational elements for 20+ bodies

### 2. **Two-Body Problem** (`twobody/`)
- **Location**: `/py-counterpart/src/poliastro/twobody/` 
- **Key Submodules**:
  - **Orbit** (`orbit/`): 1,430 lines
    - `scalar.py` (705 lines): Main `Orbit` class with 45 methods
    - `creation.py` (722 lines): 10+ `@classmethod` creation methods
  - **States** (`states.py`): 339 lines - State representation classes
  - **Elements** (`elements.py`): Core conversions
  - **Propagation** (`propagation/`): See below
  - **Angles** (`angles.py`): Anomaly conversions
  - **Sampling** (`sampling.py`): Orbit sampling utilities

### 3. **Propagation Algorithms** (`core/propagation/` + `twobody/propagation/`)
- **Location**: `/py-counterpart/src/poliastro/core/propagation/` (1,509 lines)
- **Location**: `/py-counterpart/src/poliastro/twobody/propagation/` (12 files)
- **9 Propagator Types** (with varying orbit support):
  1. **Farnocchia** (339 lines): Elliptic, parabolic, hyperbolic
  2. **Vallado** (140 lines): All orbit types
  3. **Mikkola** (131 lines): All orbit types
  4. **Pimienta** (384 lines): Elliptic, parabolic
  5. **Danby** (112 lines): All orbit types
  6. **Markley** (99 lines): Elliptic only
  7. **Gooding** (75 lines): Elliptic only
  8. **RecSeries** (122 lines): Elliptic only
  9. **Cowell** (49 lines): All orbit types + perturbations
- **Propagator Classes**: CowellPropagator, FarnocchiaPropagator, etc. with `propagate()` methods

### 4. **Core Orbital Mechanics** (`core/`)
- **Location**: `/py-counterpart/src/poliastro/core/` (7+ files)
- **Modules**:
  - `elements.py` (564 lines): 10 functions for COE ↔ RV ↔ MEE conversions
  - `angles.py` (472 lines): 15+ anomaly conversion functions
  - `iod.py` (494 lines): Initial Orbit Determination (Vallado, Izzo)
  - `maneuver.py` (manual impulses)
  - `events.py`: Event detection
  - `perturbations.py`: J2, drag, etc.

### 5. **Maneuvers** (`maneuver.py`)
- **Location**: `/py-counterpart/src/poliastro/maneuver.py` (8,033 lines)
- **Key Features**:
  - `Maneuver` class: Multiple impulses with time/velocity
  - `Hohmann`: Circular coplanar transfers
  - `Bielliptic`: For large altitude changes
  - `Lambert`: Arbitrary orbit transfers
  - `correct_pericenter`: Orbit maintenance

### 6. **Frames & Reference Systems** (`frames/`)
- **Location**: `/py-counterpart/src/poliastro/frames/` 
- **Key Components**:
  - `Planes` enum: EARTH_EQUATOR, ECLIPTIC, etc.
  - Planetary ICRSframes: 8 planet + Moon frames
  - Planetary Fixed frames: 8 planet + Moon frames
  - Ecliptic frames: GeocentricSolarEcliptic

### 7. **Plotting & Visualization** (`plotting/`)
- **Location**: `/py-counterpart/src/poliastro/plotting/` (1,500+ lines)
- **Files**:
  - `orbit/plotter.py` (807 lines): Main 2D/3D orbit visualization
  - `orbit/backends/`: Matplotlib (307), Plotly (512)
  - `gabbard.py`, `porkchop.py`, `tisserand.py`: Specialized plots
  - **Note**: Matplotlib/Plotly-dependent; lower priority for MVP

### 8. **CZML Export** (`czml/`)
- **Location**: `/py-counterpart/src/poliastro/czml/extract_czml.py` (655 lines)
- **Purpose**: Cesium JS 3D visualization format export
- **Class**: `CZMLExtractor` (45+ methods)

### 9. **Initial Orbit Determination (IOD)** (`iod/`)
- **Location**: `/py-counterpart/src/poliastro/iod/` 
- **Algorithms**:
  - Vallado IOD (iod.py in core)
  - Izzo Lambert solver (`izzo.py`)

### 10. **Ephemerides & Earth Models** (`ephem.py`, `earth/`)
- **Location**: `/py-counterpart/src/poliastro/ephem.py` (362 lines)
- **Location**: `/py-counterpart/src/poliastro/earth/` 
- **Features**:
  - Ephemerides interpolation
  - Atmosphere models (COESA 1962/1976)
  - Spheroid location calculations

### 11. **Three-Body Problem** (`threebody/`)
- **Location**: `/py-counterpart/src/poliastro/threebody/` 
- **Features**: CR3BP dynamics, Lagrange points, sphere of influence

### 12. **Math Utilities** (`_math/`)
- **Location**: `/py-counterpart/src/poliastro/_math/` 
- **Components**: linalg, integrate, ivp, optimize, special functions
- **Purpose**: Low-level numerical support

## Test Coverage Summary

| Module | Test File | Lines | Test Count |
|--------|-----------|-------|-----------|
| Orbit | `test_orbit.py` | 1,356 | 85 tests |
| Propagation | `test_propagation.py` | 476 | 15+ tests |
| Maneuvers | `test_maneuver.py` | 289 | 15+ tests |
| Events | `test_events.py` | 463 | 20+ tests |
| Angles | `test_angles.py` | 244 | 30+ tests |
| Ephemerides | `test_ephem.py` | 413 | 20+ tests |
| CZML | `test_czml.py` | 930 | Many |

---

## **PRIORITIZED MVP MIGRATION SHORTLIST FOR RUST**

### **Tier 1: Core Foundation (Must Have)**

#### **1. Bodies & Constants Module** ✓ HIGHEST PRIORITY
- **Why**: Everything depends on body definitions (Earth, Sun, etc.)
- **Scope**: ~400 lines of data + Body struct
- **Rationale**: 
  - Foundational: all orbit operations use `attractor` body
  - Low complexity: mostly data definitions with simple methods
  - Small test suite (8 tests)
- **Files to Port**:
  - `/py-counterpart/src/poliastro/bodies.py`
  - `/py-counterpart/src/poliastro/constants/general.py` (774 lines)
  - `/py-counterpart/tests/test_bodies.py`

#### **2. Core Orbital Elements & Conversions** ✓ HIGH PRIORITY
- **Why**: Foundation for all orbit representations
- **Scope**: ~560 lines core + 472 lines angles
- **Rationale**:
  - 10 core element conversion functions (COE ↔ RV ↔ MEE)
  - 15+ anomaly conversion functions (E ↔ M ↔ ν ↔ D ↔ F)
  - Well-tested with deterministic math
  - No external dependencies beyond basic math
- **Files to Port**:
  - `/py-counterpart/src/poliastro/core/elements.py`
  - `/py-counterpart/src/poliastro/core/angles.py`
  - `/py-counterpart/tests/tests_core/` (multiple test files)
- **Estimated Effort**: Low-Medium
- **Tests**: 30+ tests for angles, elements

#### **3. Orbit Class (Basic)** ✓ HIGH PRIORITY
- **Why**: Central data structure for all orbit operations
- **Scope**: 705 lines scalar.py + 722 lines creation.py
- **Rationale**:
  - Most-tested module (85 tests in test_orbit.py)
  - Wrapper around state objects with orbital properties
  - 45+ properties/methods accessing computed orbital elements
  - 10+ creation methods (from_vectors, from_classical, stationary, etc.)
- **Files to Port**:
  - `/py-counterpart/src/poliastro/twobody/orbit/scalar.py` (main Orbit class)
  - `/py-counterpart/src/poliastro/twobody/orbit/creation.py` (factory methods)
  - `/py-counterpart/src/poliastro/twobody/states.py` (state classes)
  - `/py-counterpart/tests/tests_twobody/test_orbit.py`
- **Estimated Effort**: Medium
- **Dependencies**: Bodies, Constants, Elements, Angles
- **Tests**: 85 tests

#### **4. Propagation Core (Vallado or Farnocchia)** ✓ HIGH PRIORITY
- **Why**: Essential for orbit evolution
- **Scope**: ~500 lines (one good propagator + wrapper)
- **Rationale**:
  - Start with **Farnocchia** (339 lines) or **Vallado** (140 lines)
  - Both handle elliptic, parabolic, hyperbolic orbits
  - Relatively standalone algorithms (no external solvers needed initially)
  - 15+ propagation tests
- **Files to Port**:
  - `/py-counterpart/src/poliastro/core/propagation/farnocchia.py` OR
  - `/py-counterpart/src/poliastro/core/propagation/vallado.py`
  - `/py-counterpart/src/poliastro/twobody/propagation/farnocchia.py` (wrapper)
  - `/py-counterpart/tests/tests_twobody/test_propagation.py`
- **Estimated Effort**: Medium-High (numerical methods)
- **Dependencies**: Bodies, Elements, Angles
- **Tests**: 15+ propagation tests

### **Tier 2: High-Value Add-ons (1st Phase Post-MVP)**

#### **5. Maneuvers (Hohmann, Bielliptic, Lambert)** ✓ MEDIUM-HIGH PRIORITY
- **Why**: Essential for mission design
- **Scope**: ~200 lines core logic + ~50 lines Lambert solver
- **Rationale**:
  - 15+ maneuver tests
  - Hohmann & Bielliptic are ~20 lines each
  - Lambert (Izzo) is more complex but widely used
  - High user-facing value
- **Files to Port**:
  - `/py-counterpart/src/poliastro/maneuver.py` (Maneuver class)
  - `/py-counterpart/src/poliastro/core/maneuver.py` (fast versions)
  - `/py-counterpart/src/poliastro/iod/izzo.py` (Lambert solver)
  - `/py-counterpart/tests/test_maneuver.py`
- **Estimated Effort**: Medium
- **Dependencies**: Orbit, Propagation, Elements
- **Tests**: 15+ tests

#### **6. Events (Elevation, Anomaly Crossings, Ground Track)** ✓ MEDIUM PRIORITY
- **Why**: Critical for mission analysis
- **Scope**: ~460 lines
- **Rationale**:
  - Ground station visibility crucial for ops
  - 20+ event tests
  - Integrates with propagators
- **Files to Port**:
  - `/py-counterpart/src/poliastro/core/events.py`
  - `/py-counterpart/src/poliastro/twobody/events.py`
  - `/py-counterpart/tests/tests_twobody/test_events.py`
- **Estimated Effort**: Medium
- **Dependencies**: Orbit, Propagation
- **Tests**: 20+ tests

#### **7. Frames (Basic Plane Support)** ✓ MEDIUM PRIORITY
- **Why**: Reference frame transformations needed
- **Scope**: ~250 lines (enum + basic utilities)
- **Rationale**:
  - Orbit.get_frame() needs Planes enum
  - Start with basic EARTH_EQUATOR, ECLIPTIC only
  - Planetary frames can be deferred
- **Files to Port**:
  - `/py-counterpart/src/poliastro/frames/enums.py` (Planes enum)
  - `/py-counterpart/src/poliastro/frames/util.py` (basic utilities)
  - `/py-counterpart/tests/test_frames.py`
- **Estimated Effort**: Low-Medium
- **Dependencies**: Bodies
- **Tests**: Core frame tests

### **Tier 3: Nice-to-Have / Phase 2 (Lower Priority for MVP)

#### **8. IOD (Initial Orbit Determination)** ☐ LOW-MEDIUM PRIORITY
- **Why**: Useful but not essential for basic operations
- **Scope**: 494 lines (Izzo embedded, Vallado methods)
- **Files to Port**:
  - `/py-counterpart/src/poliastro/iod/izzo.py`
  - `/py-counterpart/src/poliastro/iod/vallado.py`
  - `/py-counterpart/tests/test_iod.py`
- **Estimated Effort**: Medium-High
- **Tests**: 10+ tests
- **Defer if**: Time-constrained

#### **9. Perturbations (J2, Drag, Ephemeris)** ☐ LOWER PRIORITY
- **Why**: Needed for high-fidelity propagation
- **Scope**: ~2,500 lines (distributed)
- **Files to Port**:
  - `/py-counterpart/src/poliastro/core/perturbations.py`
  - `/py-counterpart/src/poliastro/ephem.py`
  - `/py-counterpart/src/poliastro/earth/atmosphere/`
- **Estimated Effort**: High (numerical methods + data)
- **Dependencies**: Cowell propagator (for perturbation support)
- **Defer if**: MVP focused on clean Kepler orbits

#### **10. Plotting & CZML** ☐ LOWEST PRIORITY
- **Why**: Visualization, not core mechanics
- **Scope**: 1,500+ lines but Matplotlib/Plotly dependent
- **Rationale**:
  - Depends on external visualization libraries
  - Can be replaced with user-provided plotting
  - CZML useful but not essential for mechanics
- **Defer for**: Post-MVP frontend work
- **Files**:
  - `/py-counterpart/src/poliastro/plotting/`
  - `/py-counterpart/src/poliastro/czml/`

---

## **RECOMMENDED MVP SCOPE (Tier 1 + partial Tier 2)**

### **Phase 1: Foundation (MVP Release)**
1. ✓ Bodies & Constants
2. ✓ Core Elements & Angles
3. ✓ Orbit Class
4. ✓ One Propagator (Farnocchia or Vallado)
5. ✓ Basic Maneuvers (Hohmann, Bielliptic)

**Estimated Effort**: 8-12 weeks for experienced Rust developer
**Estimated LOC**: 3,000-5,000 lines of Rust
**Test Coverage**: Ports ~200 test cases

### **Phase 2: Extended (First Update)**
6. ✓ Lambert Solver + Full Maneuvers
7. ✓ Events Detection
8. ✓ Frames (basic planes)

### **Phase 3: Advanced (Second Update)**
9. ☐ IOD
10. ☐ Perturbations + Cowell
11. ☐ Plotting/CZML (if resource permits)

---

## Key Dependencies for Rust Implementation

### **Must Have**
- Linear algebra: nalgebra, ndarray (for matrix operations)
- Special functions: Hyperbolic functions, elliptic integrals
- Root finding: For anomaly conversions, Lambert solver
- Unit handling: Custom units or uom crate
- Time: chrono crate + custom Julian date handling

### **Recommended**
- Plotting: plotly.rs, ndarray-stats
- Serialization: serde + serde_json (for CZML export)
- Performance: rayon for parallel propagation

### **Can Defer**
- Atmospheric models
- Ephemeris interpolation (user-provided initially)
- Planetary frame transformations

---

## Test File References
```
Core Tests:
  /py-counterpart/tests/tests_twobody/test_orbit.py (85 tests, 1,356 lines)
  /py-counterpart/tests/tests_twobody/test_propagation.py (15+ tests, 476 lines)
  /py-counterpart/tests/test_maneuver.py (15+ tests, 289 lines)
  /py-counterpart/tests/tests_twobody/test_angles.py (30+ tests, 244 lines)
  /py-counterpart/tests/tests_twobody/test_events.py (20+ tests, 463 lines)
  /py-counterpart/tests/test_frames.py (Core frame tests, 317 lines)
  /py-counterpart/tests/test_bodies.py (8 tests, 2,489 lines)
  /py-counterpart/tests/test_iod.py (10+ tests, 178 lines)
```

