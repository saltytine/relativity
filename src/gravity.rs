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
use crate::curvature::{metric, christoffel_symbols_spherical};
use std::f64::consts::PI;

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
        let g = metric::schwarzschild_spherical(0.0, r, theta, phi, rs);
        let gamma = christoffel_symbols_spherical(r, theta, rs);

        // geodesic equation: d^2x^mu/dtau^2 + Gamma^mu_{nu lambda} dx^nu/dtau dx^lambda/dtau = 0
        // the "force" is the negative of the connection term
        let mut force = [0.0; 4];
        let u = [vel.t, vel.x, vel.y, vel.z];
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
        }
    }
}

/// see module docs for usage and extension notes.
