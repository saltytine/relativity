# physics engine: relativistic geodesic integrator and ecs

this repository is a modular, extensible, and numerically robust 3d physics engine written in rust. it is designed for both classical and relativistic physics, with a focus on simulating geodesic motion in schwarzschild spacetime (black holes) and providing a flexible ecs (entity-component-system) foundation for future expansion.

---

## overview
this repository implements an extensible physics engine for simulating relativistic motion in curved spacetime, with a focus on schwarzschild (black hole) geodesics. the codebase is designed for both physical accuracy and numerical stability, using modern rust and an entity-component-system (ecs) architecture for flexibility and scalability.
\

## key features
- **relativistic geodesic integration:**
  - see [`src/curvature.rs`](src/curvature.rs), especially `integrate_geodesic_spherical` and `christoffel_symbols_spherical`.
  - integrates worldlines in schwarzschild spacetime using analytic christoffel symbols and a symplectic leapfrog method.
- **defensive numerics and robust integration:**
  - all coordinate singularities (e.g., r=0, θ=0/π) and floating-point nan/inf issues are handled with clamping and safe guards.
  - four-acceleration is capped to prevent runaway velocities; proper time step (`d_tau`) is adaptively reduced for large gamma (relativistic factor).
  - global gamma/energy clamping ensures all four-velocities remain physical; all normalization and construction routines are robust to NaN/Inf and unphysical states.
  - see robust coordinate conversion, normalization, and integration in [`src/curvature.rs`](src/curvature.rs) and [`src/relativity.rs`](src/relativity.rs).
- **ecs architecture:**
  - modular, scalable simulation of many particles/entities.
  - scheduler and dispatcher for running systems in order, with event system support (see [`src/ecs/scheduler/mod.rs`](src/ecs/scheduler/mod.rs)).
  - see [`src/ecs/mod.rs`](src/ecs/mod.rs) for the ecs core, and [`src/main.rs`](src/main.rs) for usage.
**extensible gravity models:**
  - Newtonian, Post-Newtonian, analytic General Relativity, and full Numerical Relativity infrastructure in [`src/gravity.rs`](src/gravity.rs).
  - See [`src/numerical_relativity.rs`](src/numerical_relativity.rs) for the grid-based metric evolution and BSSN formalism (foundation for full numerical relativity).
  - NumericalRelativityGravity now computes the metric, Christoffel symbols, and geodesic deviation using a 3D grid, with robust finite-difference derivatives and constraint damping.
  - Debugging output for metric, Christoffel symbols, and force vectors is available for diagnosing grid curvature and force calculation.
- **comprehensive four-vector algebra:**
  - lorentz boosts, minkowski and general metric inner products, four-force/acceleration, and normalization in [`src/relativity.rs`](src/relativity.rs).

---

---

## mathematical background


### 1. schwarzschild metric (spherical coordinates)
the schwarzschild metric describes spacetime around a non-rotating, uncharged black hole. the implementation is in [`metric::schwarzschild_spherical`](src/curvature.rs):

$$
  ds^2 = -\left(1 - \frac{r_s}{r}\right)c^2 dt^2 + \left(1 - \frac{r_s}{r}\right)^{-1} dr^2 + r^2 d\theta^2 + r^2 \sin^2\theta d\phi^2
$$
where $r_s = 2gm/c^2$ is the schwarzschild radius.

### 2. geodesic equation
the geodesic equation governs free-fall motion in curved spacetime. the code implements this in [`integrate_geodesic_spherical`](src/curvature.rs), using analytic christoffel symbols from [`christoffel_symbols_spherical`](src/curvature.rs):

$$
  \frac{d^2 x^\mu}{d\tau^2} + \gamma^\mu_{\nu\lambda} \frac{dx^\nu}{d\tau} \frac{dx^\lambda}{d\tau} = 0
$$
where $\gamma^\mu_{\nu\lambda}$ are the christoffel symbols (connection coefficients).

### 3. christoffel symbols (analytic)
for schwarzschild in spherical coordinates, the christoffel symbols are computed analytically in [`christoffel_symbols_spherical`](src/curvature.rs). this avoids the need for numerical derivatives and ensures physical correctness.


### 4. symplectic (leapfrog) integration
a leapfrog method is used for time evolution in [`integrate_geodesic_spherical`](src/curvature.rs), which is stable and preserves the norm of the four-velocity (timelike worldlines).

---
## implementation details

### 1. FourVector and FourVelocity
- `FourVector` (see [`src/relativity.rs`](src/relativity.rs)): represents spacetime events or four-momentum. Used throughout the codebase for all spacetime calculations.
- `FourVelocity`: Derivative of position with respect to proper time; always normalized to $-c^2$ in Minkowski metric. Construction and normalization routines are robust: all four-velocity and four-momentum creation is clamped to a global gamma/energy cap, and fallback logic ensures no unphysical (spacelike or NaN/Inf) states.
- Operations: Lorentz boosts (`lorentz_boost_x`, `lorentz_boost_beta`), inner products (`minkowski_dot`, `metric_dot`), normalization, conversion between 3-velocity and four-velocity, robust on-shell projection, and gamma clamping.

**Example:**
```rust
use relativity::{FourVector, FourVelocity};
let c = 299792458.0;
let v = FourVelocity::from_3velocity(0.1 * c, 0.0, 0.0, c);
let boosted = v.lorentz_boost_x(0.5); // boost in x-direction
```

### 2. metric and christoffel symbols
- `metric` (see [`src/curvature.rs`](src/curvature.rs)): encapsulates the metric tensor for minkowski or schwarzschild spacetime. used for all inner products and normalization.
- `christoffel_symbols_spherical`: returns analytic christoffel symbols for schwarzschild in spherical coordinates. used directly in geodesic integration.

**example:**
```rust
use curvature::{metric, christoffel_symbols_spherical};
let rs = 2.0 * g * mass / (c * c);
let gamma = christoffel_symbols_spherical(r, theta, rs);
let metric = metric::schwarzschild_spherical(t, r, theta, phi, rs);
```

### 3. geodesic and force integration
- `integrate_geodesic_spherical` (see [`src/curvature.rs`](src/curvature.rs)): Integrates a particle's worldline using the geodesic equation in spherical coordinates. Handles all coordinate singularities and normalizes the four-velocity. Used as the canonical entry point for relativistic motion.
- `FourVelocity::integrate_accel` (see [`src/relativity.rs`](src/relativity.rs)): Robustly integrates four-acceleration using RK4, with safety features:
  - Four-acceleration is capped to a maximum value to prevent runaway integration.
  - The proper time step (`d_tau`) is adaptively reduced for large gamma (relativistic factor).
  - Diagnostics are printed if capping or adaptive stepping is triggered.
  - All outputs are clamped and normalized to ensure physical validity.

**example:**
```rust
use relativity::FourVelocity;
let (x_next, u_next) = v.integrate_geodesic(&x, &metric_fn, None, mass, d_tau, dx);
let v_next = v.integrate_accel(four_accel, d_tau, c);
```

### 4. ecs (entity-component-system)
- `world` (see [`src/ecs/mod.rs`](src/ecs/mod.rs)): manages entities and their components. provides methods for entity creation, deletion, and component management.
- `physicssystem`: integrates acceleration into velocity and velocity into position for all entities. see [`physicssystem`](src/ecs/mod.rs) and its usage in [`src/main.rs`](src/main.rs).

**example:**
```rust
use ecs::{world, physicssystem};
let mut world = world::new();
world.register_component::<position>();
let entity = world.create_entity();
world.add_component(entity, position(fourvector { t: 0.0, x: 1.0, y: 0.0, z: 0.0 }));
let mut physics_system = physicssystem;
physics_system.run(&mut world);
```


### 5. gravity models

[`src/gravity.rs`](src/gravity.rs): provides several gravity models via the `GravityModel` trait:
  - `NewtonianGravity`: Classic Newtonian gravity (for reference)
  - `PostNewtonianGravity`: 1PN (weak field, slow motion) correction
  - `GeneralRelativityGravity`: Analytic Schwarzschild solution (static, spherically symmetric mass)
  - `NumericalRelativityGravity`: Full grid-based metric evolution and force calculation (see [`src/numerical_relativity.rs`](src/numerical_relativity.rs))
The `GravityKind` enum allows runtime selection between these models.

**Example:**
```rust
use gravity::{GravityKind, GravityModel, NumericalRelativityGravity};
let grid = ...; // create or load a grid
let gravity = GravityKind::NumericalRelativity(NumericalRelativityGravity::new(grid));
let force = gravity.gravity_force(&pos, &vel, m, &src_pos, &src_vel, src_m, c);
```

#### Numerical Relativity (BSSN, grid-based, with debugging)
- See [`src/numerical_relativity.rs`](src/numerical_relativity.rs) for the full infrastructure for evolving the metric on a 3D grid using the BSSN formalism.
- The `GravityKind::NumericalRelativity` variant now computes the metric, Christoffel symbols, and geodesic deviation at each grid point, and applies these to force calculations for all entities.
- The grid can be seeded with curvature (e.g., a Gaussian bump in the conformal factor), and the BSSN evolution is robustified with constraint damping and periodic boundary conditions.
- Debugging output is available: metric tensor, Christoffel symbol extrema, and force vectors are printed for each entity, allowing diagnosis of curvature and force propagation.


## numerical robustness
 all coordinate conversions and denominators are clamped to avoid division by zero (see `cartesian_to_spherical` and all uses of `.max()` and `.clamp()` in [`src/curvature.rs`](src/curvature.rs)).
 four-velocity and four-momentum construction is globally clamped to a maximum gamma/energy; normalization is checked for nan/inf and handled gracefully (see `normalize`, `clamp_gamma`, and `from_3velocity`/`from_four_momentum` in [`src/relativity.rs`](src/relativity.rs)).
 force integration (`integrate_accel`) is robust to large forces and velocities, with capping and adaptive stepping.
 collision response is robustified: post-collision four-momenta are projected on-shell and clamped, and all unphysical states are caught and corrected.
 no christoffel symbol transformation between coordinate systems (avoids known bugs and ensures physical correctness).

---
 see [`src/main.rs`](src/main.rs) for a full simulation loop. minimal example:
 ```rust
 use relativity::FourVector;
 use curvature::integrate_geodesic_spherical;
 let x = FourVector { t: 0.0, x: 10.0, y: 0.0, z: 0.0 };
 let u = FourVector { t: 1.0, x: 0.0, y: 0.1, z: 0.0 };
 let mass = 1.0;
 let g = 6.67430e-11;
 let c = 299792458.0;
 let d_tau = 0.01;
 let (x_next, u_next) = integrate_geodesic_spherical(&x, &u, mass, g, c, d_tau);
 ```

## further reading
- [carroll, s. m. (2004). *spacetime and geometry: an introduction to general relativity*](https://www.preposterousuniverse.com/spacetimeandgeometry/)
- [wald, r. m. (1984). *general relativity*](https://press.uchicago.edu/ucp/books/book/chicago/g/bo5957998.html)
- [poisson, e. (2004). *a relativist's toolkit: the mathematics of black-hole mechanics*](https://www.cambridge.org/core/books/relativists-toolkit/)
- [wikipedia: schwarzschild metric](https://en.wikipedia.org/wiki/schwarzschild_metric)
- [wikipedia: geodesic (general relativity)](https://en.wikipedia.org/wiki/geodesic_(general_relativity))

---

## how to use this repository

1. **clone and build:**
   ```sh
   git clone https://github.com/saltytine/relativity
   cd relativity
   cargo build
   ```
2. **run the example simulation:**
   ```sh
   cargo run
   ```
   this will initialize the ecs, create entities, and run a relativistic simulation loop (see [`src/main.rs`](src/main.rs)).
3. **write your own simulation:**
   - use the ecs to create entities with `position`, `velocity`, `acceleration`, `mass`, and `force` components.
   - use `integrate_geodesic_spherical` for relativistic motion, or extend the ecs with your own systems.
   - swap or extend gravity models in [`src/gravity.rs`](src/gravity.rs) as needed.
4. **testing:**
   ```sh
   cargo test
   ```
   unit tests are provided for core math and ecs features.

---

## project structure
- [`src/curvature.rs`](src/curvature.rs): metric, christoffel symbols, geodesic integrator, robust coordinate conversion.
- [`src/relativity.rs`](src/relativity.rs): four-vector algebra, four-velocity, four-force, normalization, lorentz boosts.
- [`src/gravity.rs`](src/gravity.rs): gravity models (newtonian, post-newtonian, analytic GR, and numerical relativity infrastructure), trait-based extensibility.
- [`src/numerical_relativity.rs`](src/numerical_relativity.rs): infrastructure for evolving the metric numerically on a grid (foundation for full numerical relativity).
- [`src/ecs/`](src/ecs/): ecs core (entities, components, systems, world management).
- [`src/main.rs`](src/main.rs): application entry point, ecs setup, simulation loop, and example usage.

---

## roadmap

### 1. math & physics foundation
- integrate a robust linear algebra library (`nalgebra` or `glam`).
- implement classical newtonian physics: collision detection, rigid body dynamics.
- extend to relativistic physics: time dilation, length contraction, relativistic momentum/energy, light speed constraints.
- physics integrators use 4d spacetime vectors; timestep based on proper time per object.

### 2. core engine architecture
- modular ecs (entity-component-system) foundation using a custom implementation.
- renderer abstraction layer: swap between opengl (light systems) and vulkan/metal/dx12 (heavy systems).
- clean separation of engine modules: ecs, physics, rendering, input, and tooling.

### 3. rendering
- flexible, multi-pass rendering pipeline.
- physically based rendering (pbr), shadow mapping.
- shaders for relativistic visual effects: doppler shift, aberration.
- level of detail (lod) system for performance scaling.

### 4. scaling & optimization
- multi-threading with rust async/rayon for physics and rendering.
- configurable feature toggles for low/high-end hardware.
- simd optimization and caching for hot physics loops.

### 5. tooling & api
- clean api for entity management, physics parameters, and relativity settings.
- debugging tools: spacetime interval visualization, worldlines, frame transforms.

### 6. testing & validation
- unit tests for all modules.
- physics validation against known relativistic scenarios (twin paradox, light clocks, relativistic collisions).
- benchmarking and tuning for various hardware profiles.

---


## progress tracking

### 0. numerical relativity & debugging
- [x] BSSN grid-based metric evolution (vacuum, 3+1D, robust numerics)
- [x] Christoffel symbol computation from grid metric (finite differences)
- [x] NumericalRelativityGravity: force calculation using grid Christoffel symbols
- [x] Debug output for metric, Christoffel symbols, and force at each entity
- [x] Grid seeding with curvature (Gaussian bump in conformal factor)
- [x] Constraint damping and periodic boundary conditions

### 1. core engine architecture
- [x] project structure initialized
- [x] ecs foundation (entity, component, system traits)
- [ ] ecs scheduler/dispatcher
- [x] entity management (creation, deletion)
- [x] component storage (sparse set, archetype, etc.)
- [x] system registration and execution
- [ ] event system
- [ ] module separation (ecs, physics, rendering, input, tooling)
- [ ] renderer abstraction layer (opengl, vulkan, metal, dx12)
- [ ] engine configuration system

### 2. math & physics foundation
- [ ] integrate linear algebra library (`nalgebra` or `glam`)
- [ ] vector, matrix, quaternion math utilities
- [x] classical physics engine core
  - [x] position and velocity components
  - [x] acceleration component and integration
  - [x] mass and force components (placeholders)
  - [x] physics system: acceleration → velocity → position
  - [ ] rigid body dynamics
  - [ ] collision detection (aabb, obb, sphere, mesh)
  - [ ] collision resolution
  - [x] physics timestep/integration
- [x] relativistic physics extension
  - [x] 4d spacetime vector math
  - [x] time dilation
  - [x] length contraction
  - [x] relativistic momentum & energy
  - [x] light speed constraints
  - [x] proper time-based physics integration
  - [ ] advanced physics systems (planned)

### 3. rendering
- [ ] rendering pipeline abstraction
- [ ] multi-pass rendering support
- [ ] physically based rendering (pbr)
- [ ] shadow mapping
- [ ] relativistic visual shaders (doppler, aberration)
- [ ] level of detail (lod) system
- [ ] model/mesh loader
- [ ] texture/material system
- [ ] camera system
- [ ] lighting system

### 4. scaling & optimization
- [ ] multi-threading (async, rayon)
- [ ] feature toggles/config for hardware profiles
- [ ] simd optimization for physics
- [ ] caching/reuse of calculations
- [ ] dynamic lod switching
- [ ] render resolution scaling

### 5. tooling & api
- [ ] public api for entity/component management
- [ ] physics parameter configuration
- [ ] relativity settings api
- [ ] debugging tools
  - [x] force, metric, and Christoffel symbol debugging output (NumericalRelativity)
  - [ ] spacetime interval visualization
  - [ ] worldline visualization
  - [ ] frame transform visualization
- [x] logging and diagnostics for numerical relativity

### 6. testing & validation
- [ ] unit tests for ecs
- [ ] unit tests for math utilities
- [ ] unit tests for classical physics
- [ ] unit tests for relativistic physics
- [ ] rendering validation tests
- [ ] physics scenario validation (twin paradox, light clocks, relativistic collisions)
- [ ] benchmarking on low-end hardware
- [ ] benchmarking on high-end hardware
- [ ] performance tuning and profiling

---

*this readme will be updated as the project progresses. each step will include detailed implementation, documentation, and unit tests.*
