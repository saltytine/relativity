//! # User Guide: Gravity Models and ECS Usage
//!
//! ## How to Change the Gravity Model
//!
//! 1. In your main simulation code (e.g., `main.rs`), import the gravity models:
//!    ```rust
//!    use gravity::{GravityKind, NewtonianGravity, PostNewtonianGravity, GeneralRelativityGravity, NumericalRelativityGravity};
//!    ```
//!
//! 2. Select the desired gravity model by constructing a `GravityKind`:
//!    ```rust
//!    // Newtonian gravity
//!    let gravity_model = GravityKind::Newtonian(NewtonianGravity { g_const: 6.67430e-11 });
//!    // Post-Newtonian gravity
//!    let gravity_model = GravityKind::PostNewtonian(PostNewtonianGravity { g_const: 6.67430e-11 });
//!    // Analytic General Relativity (Schwarzschild)
//!    let gravity_model = GravityKind::GeneralRelativity(GeneralRelativityGravity { g_const: 6.67430e-11 });
//!    // Numerical Relativity (infrastructure, WIP)
//!    let grid = numerical_relativity::Grid3D::new(nx, ny, nz, dx, dy, dz);
//!    let gravity_model = GravityKind::NumericalRelativity(NumericalRelativityGravity::new(grid));
//!    ```
//!
//! 3. Use `gravity_model.gravity_force(...)` wherever you compute gravitational forces between objects.
//!
//! ## How to Spawn, Collide, and Add Velocity/Acceleration to Objects
//!
//! - **Spawning Entities:**
//!   1. Use the ECS world to create a new entity:
//!      ```rust
//!      let entity = world.create_entity();
//!      ```
//!   2. Add components such as position, velocity, acceleration, and mass:
//!      ```rust
//!      world.add_component(entity, Position(FourVector { t: 0.0, x: 1.0, y: 2.0, z: 3.0 }));
//!      world.add_component(entity, Velocity(FourVelocity::from_3velocity(0.1, 0.0, 0.0, c)));
//!      world.add_component(entity, Acceleration(FourVector { t: 0.0, x: 0.0, y: 0.0, z: 0.0 }));
//!      world.add_component(entity, Mass(1.0));
//!      ```
//!
//! - **Adding/Updating Velocity and Acceleration:**
//!   - To update velocity or acceleration, use `world.add_component` with the new value. This will overwrite the previous value for that entity:
//!     ```rust
//!     world.add_component(entity, Velocity(FourVelocity::from_3velocity(0.2, 0.0, 0.0, c)));
//!     world.add_component(entity, Acceleration(FourVector { t: 0.0, x: 0.01, y: 0.0, z: 0.0 }));
//!     ```
//!
//! - **Collision Handling:**
//!   - Collisions are detected and resolved in the main simulation loop. To customize collision behavior, modify the collision detection and response code in `main.rs`.
//!   - Example: The code checks for entities within a certain distance and applies relativistic collision response, updating velocities accordingly.
//!
//! - **See Also:**
//!   - `src/main.rs` for full simulation setup and examples.
//!   - `src/ecs/mod.rs` for ECS core and system registration.
//!   - `src/relativity.rs` for four-vector and velocity math.
//!   - `src/numerical_relativity.rs` for grid-based metric evolution (numerical relativity).
//! relativistic gravity module
//!
//! this module provides a framework for advanced gravity modeling in a relativistic context
//! it supports extensible strategies: newtonian, post-newtonian, and (stub) full gr
//!
//! - for Newtonian and post-Newtonian, see the provided implementations
//! - for full gr, see the stub and documentation for future extension

use crate::relativity::{FourVector, FourVelocity, FourForce};

/// trait for gravity models
pub trait GravityModel {
    /// compute the four-force of gravity on a test mass at `pos` with velocity `vel` and mass `m`,
    /// due to a source at `src_pos` with mass `src_m` and velocity `src_vel`
    ///
    /// returns the four-force to apply to the test mass
    fn gravity_force(
        &self,
        pos: &FourVector,
        vel: &FourVelocity,
        m: f64,
        src_pos: &FourVector,
        src_vel: &FourVelocity,
        src_m: f64,
        c: f64,
    ) -> FourForce;
}

/// newtonian gravity (for reference, not relativistically correct)
pub struct NewtonianGravity {
    pub g_const: f64,
}

impl GravityModel for NewtonianGravity {
    fn gravity_force(
        &self,
        pos: &FourVector,
        _vel: &FourVelocity,
        m: f64,
        src_pos: &FourVector,
        _src_vel: &FourVelocity,
        src_m: f64,
        _c: f64,
    ) -> FourForce {
        let dx = src_pos.x - pos.x;
        let dy = src_pos.y - pos.y;
        let dz = src_pos.z - pos.z;
        let r2 = dx*dx + dy*dy + dz*dz + 1e-6;
        let r = r2.sqrt();
        let f_mag = self.g_const * m * src_m / r2;
        FourForce {
            t: 0.0,
            x: f_mag * dx / r,
            y: f_mag * dy / r,
            z: f_mag * dz / r,
        }
    }
}

/// post-Newtonian gravity (1PN, weak field, slow motion)
pub struct PostNewtonianGravity {
    pub g_const: f64,
}

impl GravityModel for PostNewtonianGravity {
    fn gravity_force(
        &self,
        pos: &FourVector,
        vel: &FourVelocity,
        m: f64,
        src_pos: &FourVector,
        src_vel: &FourVelocity,
        src_m: f64,
        c: f64,
    ) -> FourForce {
        // 1PN correction: F = F_Newton * (1 + v^2/c^2 + ...)
        let dx = src_pos.x - pos.x;
        let dy = src_pos.y - pos.y;
        let dz = src_pos.z - pos.z;
        let r2 = dx*dx + dy*dy + dz*dz + 1e-6;
        let r = r2.sqrt();
        let f_mag = self.g_const * m * src_m / r2;
        let (vx, vy, vz) = vel.three_velocity(c);
        let v2 = vx*vx + vy*vy + vz*vz;
        let pn_corr = 1.0 + v2 / (c*c); // 1PN correction (very simplified)
        FourForce {
            t: 0.0,
            x: f_mag * dx / r * pn_corr,
            y: f_mag * dy / r * pn_corr,
            z: f_mag * dz / r * pn_corr,
        }
    }
}

/// full general relativity (gr) gravity: computes the metric, christoffel symbols, and geodesic deviation
/// currently implements the schwarzschild metric for a point mass source
///
/// note: this implementation uses the analytic schwarzschild solution (static, spherically symmetric mass)
/// to handle arbitrary mass distributions, dynamic spacetimes, or to numerically solve einstein's field equations,
/// a full numerical relativity approach would be required. this is not implemented here, but the code structure
/// allows for future extension to more general metrics or numerical solutions.
/// maybe do this later
use crate::curvature::{Metric, christoffel_symbols_spherical};
use crate::numerical_relativity::NumericalRelativityGravity;


pub struct GeneralRelativityGravity {
    pub g_const: f64,
}

impl GravityModel for GeneralRelativityGravity {
    fn gravity_force(
        &self,
        pos: &FourVector,
        vel: &FourVelocity,
        m: f64,
        src_pos: &FourVector,
        _src_vel: &FourVelocity,
        src_m: f64,
        c: f64,
    ) -> FourForce {
        // interpret pos and src_pos as (t, x, y, z) in cartesian coordinates
        // convert to spherical coordinates for schwarzschild metric
        let dx = pos.x - src_pos.x;
        let dy = pos.y - src_pos.y;
        let dz = pos.z - src_pos.z;
        let r = (dx*dx + dy*dy + dz*dz).sqrt().max(1e-6);
        let theta = if r > 1e-8 { (dz / r).acos() } else { 0.0 };
        let phi = dy.atan2(dx);
        let rs = 2.0 * self.g_const * src_m / (c * c);

        // get metric and christoffel symbols at this position
        let g = Metric::schwarzschild_spherical(0.0, r, theta, phi, rs);
        let gamma = christoffel_symbols_spherical(r, theta, rs);

        // geodesic equation: d^2x^mu/dtau^2 + Gamma^mu_{nu lambda} dx^nu/dtau dx^lambda/dtau = 0
        // the "force" is the negative of the connection term
        let mut force = [0.0; 4];
        let u = [vel.0.t, vel.0.x, vel.0.y, vel.0.z];
        for mu in 0..4 {
            let mut sum = 0.0;
            for nu in 0..4 {
                for lambda in 0..4 {
                    sum += gamma[mu][nu][lambda] * u[nu] * u[lambda];
                }
            }
            force[mu] = -m * sum;
        }

        // return as FourForce (in coordinate basis)
        FourForce {
            t: force[0],
            x: force[1],
            y: force[2],
            z: force[3],
        }
    }
}

/// utility: select a gravity model
pub enum GravityKind {
    Newtonian(NewtonianGravity),
    PostNewtonian(PostNewtonianGravity),
    GeneralRelativity(GeneralRelativityGravity),
    NumericalRelativity(NumericalRelativityGravity),
}

impl GravityModel for GravityKind {
    fn gravity_force(
        &self,
        pos: &FourVector,
        vel: &FourVelocity,
        m: f64,
        src_pos: &FourVector,
        src_vel: &FourVelocity,
        src_m: f64,
        c: f64,
    ) -> FourForce {
        match self {
            GravityKind::Newtonian(g) => g.gravity_force(pos, vel, m, src_pos, src_vel, src_m, c),
            GravityKind::PostNewtonian(g) => g.gravity_force(pos, vel, m, src_pos, src_vel, src_m, c),
            GravityKind::GeneralRelativity(g) => g.gravity_force(pos, vel, m, src_pos, src_vel, src_m, c),
            GravityKind::NumericalRelativity(g) => {
                // Map the test mass position to the nearest grid point
                let (x, y, z) = (pos.x, pos.y, pos.z);
                let grid = &g.grid;
                let i = ((x / grid.dx).round() as isize).clamp(0, (grid.nx as isize) - 1) as usize;
                let j = ((y / grid.dy).round() as isize).clamp(0, (grid.ny as isize) - 1) as usize;
                let k = ((z / grid.dz).round() as isize).clamp(0, (grid.nz as isize) - 1) as usize;
                // Get the local metric and lapse/shift
                let metric = g.metric_at(i, j, k);
                let lapse = g.lapse_at(i, j, k);
                let shift = g.shift_at(i, j, k);
                // Build a 4x4 metric tensor (ADM 3+1 split)
                let mut g4 = [[0.0; 4]; 4];
                g4[0][0] = -lapse * lapse;
                for a in 0..3 {
                    g4[0][a+1] = shift[a];
                    g4[a+1][0] = shift[a];
                    for b in 0..3 {
                        g4[a+1][b+1] = metric[a][b];
                    }
                }
                // Print metric tensor for debugging
                println!("[numrel-debug] Entity at grid ({},{},{}): metric = {:?}", i, j, k, metric);
                // Compute the inverse metric
                let mut g4_inv = [[0.0; 4]; 4];
                // For simplicity, invert only the spatial part and set time-space to zero (approximate)
                let mut spatial = [[0.0; 3]; 3];
                for a in 0..3 { for b in 0..3 { spatial[a][b] = metric[a][b]; } }
                let inv_spatial = crate::numerical_relativity::invert_3x3(&spatial).unwrap_or([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
                for a in 0..3 { for b in 0..3 { g4_inv[a+1][b+1] = inv_spatial[a][b]; } }
                g4_inv[0][0] = -1.0 / (lapse * lapse + 1e-12);
                // Compute Christoffel symbols: Gamma^mu_{nu lambda} = 0.5 * g^{mu rho} (d_nu g_{rho lambda} + d_lambda g_{rho nu} - d_rho g_{nu lambda})
                let dx = grid.dx;
                let dy = grid.dy;
                let dz = grid.dz;
                let nx = grid.nx;
                let ny = grid.ny;
                let nz = grid.nz;
                // Helper to get g_{mu nu} at (i,j,k) for mu,nu in 0..4
                let get_g = |mu: usize, nu: usize, di: isize, dj: isize, dk: isize| -> f64 {
                    let ii = (i as isize + di).clamp(0, (nx as isize) - 1) as usize;
                    let jj = (j as isize + dj).clamp(0, (ny as isize) - 1) as usize;
                    let kk = (k as isize + dk).clamp(0, (nz as isize) - 1) as usize;
                    if mu == 0 && nu == 0 {
                        -g.lapse_at(ii, jj, kk).powi(2)
                    } else if mu == 0 && nu >= 1 {
                        g.shift_at(ii, jj, kk)[nu-1]
                    } else if mu >= 1 && nu == 0 {
                        g.shift_at(ii, jj, kk)[mu-1]
                    } else {
                        g.metric_at(ii, jj, kk)[mu-1][nu-1]
                    }
                };
                let mut gamma = [[[0.0; 4]; 4]; 4];
                for mu in 0..4 {
                    for nu in 0..4 {
                        for lam in 0..4 {
                            let mut sum = 0.0;
                            for rho in 0..4 {
                                // Central finite differences for derivatives
                                // For each derivative, select direction and spacing
                                let (d_nu, h_nu) = match nu {
                                    1 => ((1,0,0), dx),
                                    2 => ((0,1,0), dy),
                                    3 => ((0,0,1), dz),
                                    _ => ((0,0,0), 1.0),
                                };
                                let (d_lam, h_lam) = match lam {
                                    1 => ((1,0,0), dx),
                                    2 => ((0,1,0), dy),
                                    3 => ((0,0,1), dz),
                                    _ => ((0,0,0), 1.0),
                                };
                                let (d_rho, h_rho) = match rho {
                                    1 => ((1,0,0), dx),
                                    2 => ((0,1,0), dy),
                                    3 => ((0,0,1), dz),
                                    _ => ((0,0,0), 1.0),
                                };
                                let dg_rho_lam_dnu = (get_g(rho, lam, d_nu.0, d_nu.1, d_nu.2) - get_g(rho, lam, -(d_nu.0), -(d_nu.1), -(d_nu.2))) / (2.0 * h_nu);
                                let dg_rho_nu_dlam = (get_g(rho, nu, d_lam.0, d_lam.1, d_lam.2) - get_g(rho, nu, -(d_lam.0), -(d_lam.1), -(d_lam.2))) / (2.0 * h_lam);
                                let dg_nu_lam_drho = (get_g(nu, lam, d_rho.0, d_rho.1, d_rho.2) - get_g(nu, lam, -(d_rho.0), -(d_rho.1), -(d_rho.2))) / (2.0 * h_rho);
                                sum += g4_inv[mu][rho] * (dg_rho_lam_dnu + dg_rho_nu_dlam - dg_nu_lam_drho);
                            }
                            gamma[mu][nu][lam] = 0.5 * sum;
                        }
                    }
                }
                // Print a summary of Christoffel symbols for debugging
                let mut max_gamma = 0.0;
                let mut min_gamma = 0.0;
                for mu in 0..4 { for nu in 0..4 { for lam in 0..4 {
                    let val = gamma[mu][nu][lam];
                    if val > max_gamma { max_gamma = val; }
                    if val < min_gamma { min_gamma = val; }
                }}}
                println!("[numrel-debug] Christoffel symbols at ({},{},{}): max = {:.3e}, min = {:.3e}", i, j, k, max_gamma, min_gamma);
                // Use the four-velocity of the test mass
                let u = [vel.0.t, vel.0.x, vel.0.y, vel.0.z];
                let mut force = [0.0; 4];
                for mu in 0..4 {
                    let mut sum = 0.0;
                    for nu in 0..4 {
                        for lambda in 0..4 {
                            sum += gamma[mu][nu][lambda] * u[nu] * u[lambda];
                        }
                    }
                    force[mu] = -m * sum;
                }
                println!("[numrel-debug] FourForce at ({},{},{}): [{:.3e}, {:.3e}, {:.3e}, {:.3e}]", i, j, k, force[0], force[1], force[2], force[3]);
                FourForce {
                    t: force[0],
                    x: force[1],
                    y: force[2],
                    z: force[3],
                }
            }
        }
    }
}

