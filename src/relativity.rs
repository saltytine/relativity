impl FourVector {
    /// clamp the spatial part so that gamma (t/c) does not exceed max_gamma
    pub fn clamp_gamma(&self, c: f64, max_gamma: f64) -> FourVector {
        let gamma = self.t.abs() / c;
        if gamma > max_gamma {
            let spatial2 = self.x * self.x + self.y * self.y + self.z * self.z;
            let max_spatial2 = (max_gamma * max_gamma - 1.0) * c * c;
            let scale = if spatial2 > 0.0 && spatial2 > max_spatial2 {
                (max_spatial2 / spatial2).sqrt()
            } else {
                1.0
            };
            eprintln!(
                "[warn] FourVector::clamp_gamma: gamma = {:.3e} > max_gamma = {:.3e}, clamping spatial part.\n  Before: t = {:.3e}, |p| = {:.3e}\n  After:  t = {:.3e}, |p| = {:.3e}",
                gamma, max_gamma, self.t, spatial2.sqrt(), self.t, (spatial2 * scale * scale).sqrt()
            );
            FourVector {
                t: self.t,
                x: self.x * scale,
                y: self.y * scale,
                z: self.z * scale,
            }
        } else {
            *self
        }
    }
    /// project this four-momentum onto the mass shell: set t = sqrt(|p|^2 + (m c)^2)
    /// if t is negative, keeps the sign of the original t (for future-proofing, e.g. antiparticles)
    pub fn onshell_time(&self, mass: f64, c: f64) -> f64 {
        let spatial2 = self.x * self.x + self.y * self.y + self.z * self.z;
        let mc2 = mass * c;
        let t_sign = if self.t >= 0.0 { 1.0 } else { -1.0 };
        let t_sq = spatial2 + mc2 * mc2;
        if t_sq < 0.0 {
            eprintln!(
                "[fatal] FourVector::onshell_time: negative argument to sqrt! spatial2 = {:.6}, mc2^2 = {:.6}, t_sq = {:.6}",
                spatial2, mc2 * mc2, t_sq
            );
            0.0
        } else {
            t_sign * t_sq.sqrt()
        }
    }

    /// return a new FourVector with the same spatial part, but t set to the on-shell value
    pub fn onshell_project(&self, mass: f64, c: f64) -> FourVector {
        let t_new = self.onshell_time(mass, c);
        FourVector { t: t_new, x: self.x, y: self.y, z: self.z }
    }
}
use std::ops::{Sub, Add};
use crate::curvature::Metric;
// implement subtraction for FourVector
impl Sub for FourVector {
    type Output = FourVector;
    fn sub(self, rhs: FourVector) -> FourVector {
        FourVector {
            t: self.t - rhs.t,
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

// implement addition for FourVector (optional, for symmetry)
impl Add for FourVector {
    type Output = FourVector;
    fn add(self, rhs: FourVector) -> FourVector {
        FourVector {
            t: self.t + rhs.t,
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}
/// four-vector for spacetime (ct, x, y, z)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FourVector {
    pub t: f64, // ct
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl FourVector {
    /// general Lorentz boost by a 3D beta vector (beta = v / c).
    /// beta components must satisfy beta² < 1 (we check and clamp).
    pub fn lorentz_boost_beta(&self, beta_x: f64, beta_y: f64, beta_z: f64, c: f64) -> FourVector {
        let beta2 = beta_x*beta_x + beta_y*beta_y + beta_z*beta_z;
        if beta2 <= 0.0 {
            debug_assert!(self.t.is_finite() && self.x.is_finite() && self.y.is_finite() && self.z.is_finite(), "FourVector contains NaN or Inf");
            return *self;
        }
        // safety: clamp beta2 slightly below 1
        let eps = 1e-12;
        let max_beta2 = 1.0 - eps;
        let beta2_clamped = if beta2 >= max_beta2 { max_beta2 } else { beta2 };
        debug_assert!(beta2_clamped >= 0.0 && beta2_clamped < 1.0, "beta2_clamped out of range");
        let gamma = 1.0 / (1.0 - beta2_clamped).sqrt();

        // spatial r vector
        let r_x = self.x;
        let r_y = self.y;
        let r_z = self.z;

        // beta·r
        let beta_dot_r = beta_x * r_x + beta_y * r_y + beta_z * r_z;
        // t is ct (length units)
        let t_ct = self.t;

        // transformed time (ct')
        let t_prime = gamma * (t_ct - beta_dot_r / c);

        // compute factor used for spatial transform:
        // r' = r + ( (γ - 1) * (β·r)/β² - γ * t ) * β
        let factor = if beta2_clamped > 0.0 {
            ((gamma - 1.0) * beta_dot_r / beta2_clamped) - gamma * t_ct
        } else {
            -gamma * t_ct
        };

        let rpx = r_x + beta_x * factor;
        let rpy = r_y + beta_y * factor;
        let rpz = r_z + beta_z * factor;

        let out = FourVector {
            t: t_prime,
            x: rpx,
            y: rpy,
            z: rpz,
        };
        debug_assert!(out.t.is_finite() && out.x.is_finite() && out.y.is_finite() && out.z.is_finite(), "FourVector contains NaN or Inf after boost");
        out
    }
    /// minkowski inner product (metric signature: -+++)
    pub fn minkowski_dot(self, other: FourVector) -> f64 {
        -self.t * other.t + self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// general inner product with arbitrary metric (default: Minkowski)
    pub fn metric_dot(self, other: FourVector, metric: &Metric) -> f64 {
        let a = [self.t, self.x, self.y, self.z];
        let b = [other.t, other.x, other.y, other.z];
        let mut sum = 0.0;
        for mu in 0..4 {
            for nu in 0..4 {
                sum += metric.g[mu][nu] * a[mu] * b[nu];
            }
        }
        sum
    }

    /// spacetime interval squared (s^2 = -c^2 t^2 + x^2 + y^2 + z^2)
    pub fn interval2(self) -> f64 {
        self.minkowski_dot(self)
    }

    /// spacetime interval squared with arbitrary metric
    pub fn interval2_metric(self, metric: &Metric) -> f64 {
        self.metric_dot(self, metric)
    }

    /// lorentz boost in x direction (beta = v/c)
    pub fn lorentz_boost_x(self, beta: f64) -> FourVector {
        let gamma = 1.0 / (1.0 - beta * beta).sqrt();
        FourVector {
            t: gamma * (self.t - beta * self.x),
            x: gamma * (self.x - beta * self.t),
            y: self.y,
            z: self.z,
        }
    }

    /// add two four-vectors
    pub fn add(self, other: FourVector) -> FourVector {
        FourVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    /// scale a four-vector
    pub fn scale(self, s: f64) -> FourVector {
        FourVector {
            t: self.t * s,
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
}

/// 4-velocity: derivative of position four-vector wrt proper time

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FourVelocity(pub FourVector);

impl FourVelocity {
    /// general Lorentz boost by a 3D beta vector (beta = v / c)
    pub fn lorentz_boost_beta(&self, bx: f64, by: f64, bz: f64, c: f64) -> FourVelocity {
        FourVelocity(self.0.lorentz_boost_beta(bx, by, bz, c))
    }
    /// compute the four-momentum for a given mass and c
    pub fn four_momentum(&self, mass: f64, c: f64) -> FourVector {
        FourVector {
            t: self.0.t * mass,
            x: self.0.x * mass,
            y: self.0.y * mass,
            z: self.0.z * mass,
        }
    }

    /// construct a FourVelocity from a four-momentum, mass, and c
    pub fn from_four_momentum(p: FourVector, mass: f64, c: f64) -> FourVelocity {
        // always project the four-momentum onto the mass shell before constructing the four-velocity
        let max_gamma = 1e4;
        let fv_proj = FourVector {
            t: p.t / mass,
            x: p.x / mass,
            y: p.y / mass,
            z: p.z / mass,
        };
        let fv_clamped = fv_proj.clamp_gamma(c, max_gamma);
        let fv_onshell = fv_clamped.onshell_project(1.0, c); // fv_proj is already divided by mass, so use m=1
        let norm2 = fv_onshell.minkowski_dot(fv_onshell);
        let tol = (c * c) * 1e-12;
        if norm2 >= -tol {
            eprintln!(
                "[fatal] FourVelocity::from_four_momentum: on-shell projection still spacelike or lightlike!\n  p = {:?}, mass = {}, c = {}\n  fv_proj = {:?}, fv_clamped = {:?}, fv_onshell = {:?}, norm2 = {:.6}, tol = {:.6}",
                p, mass, c, fv_proj, fv_clamped, fv_onshell, norm2, tol
            );
            // fallback: at-rest four-velocity
            FourVelocity(FourVector { t: c, x: 0.0, y: 0.0, z: 0.0 })
        } else {
            FourVelocity(fv_onshell)
        }
    }
    /// construct from 3-velocity (vx, vy, vz) in absolute units (same units as c)
    /// vx, vy, vz, and c must all be in the same units (e.g., m/s). Not fractions of c.
    pub fn from_3velocity(vx: f64, vy: f64, vz: f64, c: f64) -> FourVelocity {
        let v2 = vx * vx + vy * vy + vz * vz;
        let v2 = v2.min(c * c * 0.999_999); // clamp to just below c
        let gamma = 1.0 / (1.0 - v2 / (c * c)).sqrt();
        let max_gamma = 1e4;
        let gamma_clamped = gamma.min(max_gamma);
        if gamma > max_gamma {
            let scale = (1.0 - 1.0 / (max_gamma * max_gamma)).sqrt() * c / (v2.sqrt().max(1e-12));
            eprintln!(
                "[warn] FourVelocity::from_3velocity: gamma = {:.3e} > max_gamma = {:.3e}, clamping velocity.\n  Before: |v| = {:.3e}, After: |v| = {:.3e}",
                gamma, max_gamma, v2.sqrt(), scale * v2.sqrt()
            );
            FourVelocity(FourVector {
                t: gamma_clamped * c,
                x: vx * scale,
                y: vy * scale,
                z: vz * scale,
            })
        } else {
            FourVelocity(FourVector {
                t: gamma * c,
                x: gamma * vx,
                y: gamma * vy,
                z: gamma * vz,
            })
        }
    }

    /// get the 3-velocity (vx, vy, vz) from the four-velocity
    pub fn three_velocity(&self, c: f64) -> (f64, f64, f64) {
        let gamma = self.0.t / c;
        if gamma > 0.0 {
            (self.0.x / gamma, self.0.y / gamma, self.0.z / gamma)
        } else {
            (0.0, 0.0, 0.0)
        }
    }

    /// gamma factor (Lorentz factor)
    pub fn gamma(&self, c: f64) -> f64 {
        self.0.t / c
    }

    /// minkowski norm squared (should be -c^2 for 4-velocity)
    pub fn norm2(&self, c: f64) -> f64 {
        self.0.minkowski_dot(self.0)
    }

    /// normalize the four-velocity to -c^2 Minkowski norm
    /// panics if norm is too close to zero (indicates a bug or unphysical state)
    pub fn normalize(&self, c: f64) -> FourVelocity {
        let norm = self.norm2(c);
        let tol = (c * c) * 1e-12;
        if norm < -tol {
            let scale = (-c * c / norm).sqrt();
            FourVelocity(FourVector {
                t: self.0.t * scale,
                x: self.0.x * scale,
                y: self.0.y * scale,
                z: self.0.z * scale,
            })
        } else {
            eprintln!(
                "[fatal] FourVelocity::normalize: norm not negative and large enough!\n  self = {:?}, norm = {:.12}, tol = {:.12}",
                self, norm, tol
            );
            // fallback: at-rest four-velocity
            FourVelocity(FourVector { t: c, x: 0.0, y: 0.0, z: 0.0 })
        }
    }

    /// lorentz boost in x direction (beta = v/c)
    pub fn lorentz_boost_x(&self, beta: f64) -> FourVelocity {
        FourVelocity(self.0.lorentz_boost_x(beta))
    }

    /// general geodesic integrator (Schwarzschild, spherical coordinates)
    pub fn integrate_geodesic(
        &self,
        x: &FourVector,
        _metric_at: &dyn Fn(&FourVector) -> Metric, // unused
        _force: Option<&FourVector>, // unused
        mass: f64,
        d_tau: f64,
        _dx: f64, // unused
    ) -> (FourVector, FourVelocity) {
        let G = 6.67430e-11; // gravitational constant in m³/kg/s²
        let c = 299792458.0; // speed of light in m/s
        let (x_next, u_next) = crate::curvature::integrate_geodesic_spherical(x, &self.0, mass, G, c, d_tau);
        (x_next, FourVelocity(u_next))
    }

    /// Integrate four-acceleration using RK4 for stability
    pub fn integrate_accel(&self, four_accel: FourVector, d_tau: f64, c: f64) -> FourVelocity {
        // helper: returns four-acceleration at (u, a)
        fn four_accel_fn(u: &FourVector, a: &FourVector, c: f64) -> FourVector {
            // project a orthogonal to u (enforce u·a = 0)
            let dot = -u.t * a.t + u.x * a.x + u.y * a.y + u.z * a.z;
            let norm = -u.t * u.t + u.x * u.x + u.y * u.y + u.z * u.z;
            let factor = if norm.abs() > 1e-6 { dot / norm } else { 0.0 };
            FourVector {
                t: a.t - factor * u.t,
                x: a.x - factor * u.x,
                y: a.y - factor * u.y,
                z: a.z - factor * u.z,
            }
        }

        // cap four-acceleration magnitude to avoid runaway integration
        let mut capped_accel = four_accel;
        let accel_spatial2 = four_accel.x * four_accel.x + four_accel.y * four_accel.y + four_accel.z * four_accel.z;
        let accel_max = 10.0 * c; // arbitrary: 10x speed of light per unit proper time
        let accel_norm = accel_spatial2.sqrt();
        if accel_norm > accel_max {
            let scale = accel_max / accel_norm;
            eprintln!(
                "[safety] Capping four-acceleration: |a| = {:.3e} > accel_max = {:.3e}, scaling spatial part.",
                accel_norm, accel_max
            );
            capped_accel.x *= scale;
            capped_accel.y *= scale;
            capped_accel.z *= scale;
        }

        // adaptively reduce d_tau for large velocities
        let gamma = self.0.t.abs() / c;
        let mut d_tau_eff = d_tau;
        if gamma > 100.0 {
            d_tau_eff = d_tau / (gamma / 100.0);
            eprintln!(
                "[safety] Reducing d_tau for large gamma: gamma = {:.3e}, d_tau = {:.3e} -> {:.3e}",
                gamma, d_tau, d_tau_eff
            );
        }

        let u0 = self.0;
        let a0 = four_accel_fn(&u0, &capped_accel, c);
        let u1 = FourVector {
            t: u0.t + 0.5 * d_tau_eff * a0.t,
            x: u0.x + 0.5 * d_tau_eff * a0.x,
            y: u0.y + 0.5 * d_tau_eff * a0.y,
            z: u0.z + 0.5 * d_tau_eff * a0.z,
        };
        let a1 = four_accel_fn(&u1, &capped_accel, c);
        let u2 = FourVector {
            t: u0.t + 0.5 * d_tau_eff * a1.t,
            x: u0.x + 0.5 * d_tau_eff * a1.x,
            y: u0.y + 0.5 * d_tau_eff * a1.y,
            z: u0.z + 0.5 * d_tau_eff * a1.z,
        };
        let a2 = four_accel_fn(&u2, &capped_accel, c);
        let u3 = FourVector {
            t: u0.t + d_tau_eff * a2.t,
            x: u0.x + d_tau_eff * a2.x,
            y: u0.y + d_tau_eff * a2.y,
            z: u0.z + d_tau_eff * a2.z,
        };
        let a3 = four_accel_fn(&u3, &capped_accel, c);

        // RK4 step
        let u_next = FourVector {
            t: u0.t + (d_tau_eff / 6.0) * (a0.t + 2.0 * a1.t + 2.0 * a2.t + a3.t),
            x: u0.x + (d_tau_eff / 6.0) * (a0.x + 2.0 * a1.x + 2.0 * a2.x + a3.x),
            y: u0.y + (d_tau_eff / 6.0) * (a0.y + 2.0 * a1.y + 2.0 * a2.y + a3.y),
            z: u0.z + (d_tau_eff / 6.0) * (a0.z + 2.0 * a1.z + 2.0 * a2.z + a3.z),
        };
        // clamp gamma before normalization
        let max_gamma = 1e4;
        let u_next_clamped = u_next.clamp_gamma(c, max_gamma);
        // renormalize to -c^2 minkowski norm
        let norm = u_next_clamped.minkowski_dot(u_next_clamped);
        let tol = (c * c) * 1e-12;
        if norm < -tol {
            let scale = (-c * c / norm).sqrt();
            FourVelocity(FourVector {
                t: u_next_clamped.t * scale,
                x: u_next_clamped.x * scale,
                y: u_next_clamped.y * scale,
                z: u_next_clamped.z * scale,
            })
        } else {
            eprintln!(
                "[fatal] FourVelocity::integrate_accel: norm not negative and large enough!\n  self = {:?}, four_accel = {:?}, d_tau = {:.6}, c = {:.6}, u_next = {:?}, norm = {:.12}, tol = {:.12}",
                self, four_accel, d_tau_eff, c, u_next_clamped, norm, tol
            );
            // fallback: at-rest four-velocity
            FourVelocity(FourVector { t: c, x: 0.0, y: 0.0, z: 0.0 })
        }
    }
}

/// four-force: change in four-momentum per unit proper time
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FourForce {
    pub t: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl FourForce {
    /// project four-force orthogonal to four-velocity (SR gravity, EM, etc)
    pub fn orthogonal_to(self, u: FourVelocity, c: f64) -> FourForce {
        // F^μ = F^μ - (F_ν u^ν) u^μ / c^2
        let dot = -self.t * u.0.t + self.x * u.0.x + self.y * u.0.y + self.z * u.0.z;
        let factor = dot / (c * c);
        FourForce {
            t: self.t - factor * u.0.t,
            x: self.x - factor * u.0.x,
            y: self.y - factor * u.0.y,
            z: self.z - factor * u.0.z,
        }
    }
}
