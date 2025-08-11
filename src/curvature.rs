macro_rules! assert_finite {
    ($val:expr, $msg:expr) => {
        if !$val.is_finite() {
            eprintln!("[NaN/Inf DETECTED] {}: value = {:?}", $msg, $val);
            panic!("[NaN/Inf DETECTED] {}: value = {:?}", $msg, $val);
        }
    };
}

use crate::relativity::FourVector;

/// a 4x4 symmetric metric tensor g_{μν} at a spacetime point
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Metric {
    pub g: [[f64; 4]; 4],
}

impl Metric {
    /// minkowski metric (flat spacetime)
    pub fn minkowski() -> Self {
        Metric {
            g: [
                [-1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// schwarzschild metric in spherical (t, r, θ, φ)
    pub fn schwarzschild_spherical(_t: f64, r: f64, theta: f64, _phi: f64, rs: f64) -> Self {
        // note: i didnt try to be clever at r <= rs here, caller should guard
        let f = (1.0 - rs / r).max(1e-16);
        Metric {
            g: [
                [-f, 0.0, 0.0, 0.0],
                [0.0, 1.0 / f, 0.0, 0.0],
                [0.0, 0.0, r * r, 0.0],
                [0.0, 0.0, 0.0, r * r * theta.sin() * theta.sin()],
            ],
        }
    }
}

/// compute christoffel symbols in spherical coordinates (schwarzschild)
pub fn christoffel_symbols_spherical(r: f64, theta: f64, rs: f64) -> [[[f64; 4]; 4]; 4] {
    let mut gamma = [[[0.0; 4]; 4]; 4];
    let eps = 1e-16;
    if r <= 0.0 {
        return gamma;
    }
    let f = 1.0 - rs / r;
    let f_prime = rs / (r * r);
    let sin_theta = theta.sin();
    let cos_theta = theta.cos();
    let near_pole = sin_theta.abs() < eps;

    // time components
    if f.abs() > eps {
        gamma[0][0][1] = f_prime / (2.0 * f); // Γ^t_{tr} = Γ^t_{rt}
        gamma[0][1][0] = gamma[0][0][1];
    }

    // radial components
    if f.abs() > eps {
        gamma[1][0][0] = f * f_prime / 2.0; // Γ^r_{tt}
        gamma[1][1][1] = -f_prime / (2.0 * f); // Γ^r_{rr}
    }
    gamma[1][2][2] = -r * f; // Γ^r_{θθ}
    gamma[1][3][3] = -r * f * sin_theta * sin_theta; // Γ^r_{φφ}

    // angular components
    gamma[2][1][2] = 1.0 / r; // Γ^θ_{rθ} = Γ^θ_{θr}
    gamma[2][2][1] = gamma[2][1][2];
    gamma[2][3][3] = -sin_theta * cos_theta; // Γ^θ_{φφ}

    gamma[3][1][3] = 1.0 / r; // Γ^φ_{rφ} = Γ^φ_{φr}
    gamma[3][3][1] = gamma[3][1][3];
    // cotθ = cosθ / sinθ, but set to 0 at poles
    gamma[3][2][3] = if near_pole { 0.0 } else { cos_theta / sin_theta };
    gamma[3][3][2] = gamma[3][2][3];
    gamma
}

/// solve 3x3 linear system A x = b, returns Option<[x,y,z]>
fn solve_3x3(a: [[f64; 3]; 3], b: [f64; 3], tol: f64) -> Option<[f64; 3]> {
    // determinant
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

    if det.abs() <= tol {
        return None;
    }
    let inv_det = 1.0 / det;

    // adjugate / inverse * b
    let inv = [
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det,
            -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) * inv_det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det,
        ],
        [
            -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) * inv_det,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det,
            -(a[0][0] * a[1][2] - a[0][2] * a[1][0]) * inv_det,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det,
            -(a[0][0] * a[2][1] - a[0][1] * a[2][0]) * inv_det,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det,
        ],
    ];

    let x = [
        inv[0][0] * b[0] + inv[0][1] * b[1] + inv[0][2] * b[2],
        inv[1][0] * b[0] + inv[1][1] * b[1] + inv[1][2] * b[2],
        inv[2][0] * b[0] + inv[2][1] * b[1] + inv[2][2] * b[2],
    ];
    Some(x)
}

/// integrate geodesic in schwarzschild spacetime using spherical coordinates
/// x and u are given in Cartesian basis (FourVector::t,x,y,z).
pub fn integrate_geodesic_spherical(
    x: &FourVector, // cartesian input
    u: &FourVector, // cartesian input
    mass: f64,
    G: f64,
    c: f64,
    d_tau: f64,
) -> (FourVector, FourVector) {
    let tol = 1e-12;
    let rs = 2.0 * G * mass / (c * c);


    // convert position to spherical
    let (t0, mut r, mut theta, mut phi) = Metric::cartesian_to_spherical(x);
    assert_finite!(t0, "t0 (cartesian_to_spherical)");
    assert_finite!(r, "r (cartesian_to_spherical)");
    assert_finite!(theta, "theta (cartesian_to_spherical)");
    assert_finite!(phi, "phi (cartesian_to_spherical)");
    if r <= rs + tol {
        // refuse to step across horizon in these coords
        return (*x, *u);
    }


    // precompute trig
    let st = theta.sin();
    let ct = theta.cos();
    let sp = phi.sin();
    let cp = phi.cos();
    let near_pole = st.abs() < tol;
    assert_finite!(st, "sin(theta)");
    assert_finite!(ct, "cos(theta)");
    assert_finite!(sp, "sin(phi)");
    assert_finite!(cp, "cos(phi)");

    // build Jacobian matrix A (3x3): ∂(x,y,z)/∂(r,θ,φ) columns
    // A * [dr, dθ, dφ]^T = [u_x, u_y, u_z]^T
    let a = [
        [st * cp,  r * ct * cp,   -r * st * sp],
        [st * sp,  r * ct * sp,    r * st * cp],
        [ct,       -r * st,        0.0         ],
    ];

    assert_finite!(u.x, "u.x (input)");
    assert_finite!(u.y, "u.y (input)");
    assert_finite!(u.z, "u.z (input)");
    let b = [u.x, u.y, u.z];

    // solve for dr, dθ, dφ
    let solved = solve_3x3(a, b, 1e-20);
    let (dr_dtau, dtheta_dtau, dphi_dtau) = match solved {
        Some(sol) => {
            assert_finite!(sol[0], "dr_dtau (solve_3x3)");
            assert_finite!(sol[1], "dtheta_dtau (solve_3x3)");
            let dphi = if near_pole { 0.0 } else { sol[2] };
            assert_finite!(dphi, "dphi_dtau (solve_3x3, pole-guarded)");
            (sol[0], sol[1], dphi)
        },
        None => {
            // degenerate geometry (e.g., st ~ 0 or r ~ 0)
            // fall back to guarded special-cases.
            // if at pole, set dphi = 0, cotθ = 0
            if near_pole {
                // reduced 2x2: columns for r and θ (ignore φ)
                let a2 = [
                    [st * cp, r * ct * cp],
                    [st * sp, r * ct * sp],
                ];
                let b2 = [u.x, u.y];
                let det2 = a2[0][0] * a2[1][1] - a2[0][1] * a2[1][0];
                if det2.abs() < 1e-20 {
                    // disaster: fall back to radial-only approximation
                    let dr = (x.x * u.x + x.y * u.y + x.z * u.z) / r;
                    assert_finite!(dr, "dr (radial fallback)");
                    (dr, 0.0, 0.0)
                } else {
                    let inv_det2 = 1.0 / det2;
                    let dr = inv_det2 * (b2[0] * a2[1][1] - a2[0][1] * b2[1]);
                    let dtheta = inv_det2 * (a2[0][0] * b2[1] - b2[0] * a2[1][0]);
                    assert_finite!(dr, "dr (2x2 fallback)");
                    assert_finite!(dtheta, "dtheta (2x2 fallback)");
                    (dr, dtheta, 0.0)
                }
            } else {
                // fallback: radial projection only
                let dr = (x.x * u.x + x.y * u.y + x.z * u.z) / r;
                assert_finite!(dr, "dr (radial fallback 2)");
                (dr, 0.0, 0.0)
            }
        }
    };


    // build u_sph (coordinate basis rates)
    assert_finite!(u.t, "u.t (input)");
    let u_sph = [u.t, dr_dtau, dtheta_dtau, dphi_dtau];


    // christoffel symbols in spherical
    let f = 1.0 - rs / r;
    if f.abs() < tol {
        // schwarzschild f near 0: horizon, step size must shrink or bail
        eprintln!("[geodesic] Schwarzschild f near 0 (r ≈ r_s): integration step too large or at horizon. r = {}, r_s = {}", r, rs);
        return (*x, *u);
    }
    let gamma = christoffel_symbols_spherical(r, theta, rs);
    for mu in 0..4 { for nu in 0..4 { for lam in 0..4 {
        assert_finite!(gamma[mu][nu][lam], "Christoffel symbol");
    }}}


    // leapfrog-style integration in spherical coords
    let mut x_sph = [t0, r, theta, phi];
    let mut u_sph_next = u_sph;
    for i in 0..4 { assert_finite!(x_sph[i], "x_sph (init)"); assert_finite!(u_sph[i], "u_sph (init)"); }


    // half step in position
    for i in 0..4 {
        x_sph[i] += 0.5 * d_tau * u_sph[i];
        assert_finite!(x_sph[i], "x_sph (half step)");
    }


    // full step in velocity: u_next = u + dτ * (-Γ^μ_{νλ} u^ν u^λ)
    for mu in 0..4 {
        let mut acc = 0.0;
        for nu in 0..4 {
            for lam in 0..4 {
                acc -= gamma[mu][nu][lam] * u_sph[nu] * u_sph[lam];
            }
        }
        assert_finite!(acc, "geodesic acc");
        u_sph_next[mu] += d_tau * acc;
        assert_finite!(u_sph_next[mu], "u_sph_next (after acc)");
    }


    // normalize four-velocity in spherical metric
    let g_sph = Metric::schwarzschild_spherical(x_sph[0], x_sph[1], x_sph[2], x_sph[3], rs);
    let mut norm2 = 0.0;
    let u_arr = [u_sph_next[0], u_sph_next[1], u_sph_next[2], u_sph_next[3]];
    for mu in 0..4 {
        for nu in 0..4 {
            norm2 += g_sph.g[mu][nu] * u_arr[mu] * u_arr[nu];
        }
    }
    assert_finite!(norm2, "norm2 (before normalization)");

    // only normalize if clearly timelike (norm2 < -eps_norm)
    let eps_norm = 1e-12;
    if norm2 < -eps_norm {
        let scale = (-1.0 / norm2).sqrt();
        assert_finite!(scale, "scale (normalization)");
        for i in 0..4 {
            u_sph_next[i] *= scale;
            assert_finite!(u_sph_next[i], "u_sph_next (after normalization)");
        }
    } // do not abs() norm2, do not normalize if not timelike


    // final half-step in position
    for i in 0..4 {
        x_sph[i] += 0.5 * d_tau * u_sph_next[i];
        assert_finite!(x_sph[i], "x_sph (final half step)");
    }


    // if new radius is at/inside horizon, clamp and return previous
    let r_new = x_sph[1];
    assert_finite!(r_new, "r_new (after integration)");
    if r_new <= rs + 1e-10 {
        return (*x, *u);
    }

    // convert back to Cartesian coordinates
    let t_new = x_sph[0];
    let r_new = x_sph[1];
    let theta_new = x_sph[2];
    let phi_new = x_sph[3];


    let st = theta_new.sin();
    let ct = theta_new.cos();
    let sp = phi_new.sin();
    let cp = phi_new.cos();
    assert_finite!(st, "sin(theta_new)");
    assert_finite!(ct, "cos(theta_new)");
    assert_finite!(sp, "sin(phi_new)");
    assert_finite!(cp, "cos(phi_new)");


    let x_cart = FourVector {
        t: t_new,
        x: r_new * st * cp,
        y: r_new * st * sp,
        z: r_new * ct,
    };
    assert_finite!(x_cart.t, "x_cart.t");
    assert_finite!(x_cart.x, "x_cart.x");
    assert_finite!(x_cart.y, "x_cart.y");
    assert_finite!(x_cart.z, "x_cart.z");

    // convert velocities: dx/dτ = A * [dr, dθ, dφ]
    let a_new = [
        [st * cp,  r_new * ct * cp,   -r_new * st * sp],
        [st * sp,  r_new * ct * sp,    r_new * st * cp],
        [ct,       -r_new * st,        0.0         ],
    ];
    let dr = u_sph_next[1];
    let dth = u_sph_next[2];
    let dph = u_sph_next[3];


    let vx = a_new[0][0] * dr + a_new[0][1] * dth + a_new[0][2] * dph;
    let vy = a_new[1][0] * dr + a_new[1][1] * dth + a_new[1][2] * dph;
    let vz = a_new[2][0] * dr + a_new[2][1] * dth + a_new[2][2] * dph;
    assert_finite!(vx, "vx (cartesian velocity)");
    assert_finite!(vy, "vy (cartesian velocity)");
    assert_finite!(vz, "vz (cartesian velocity)");

    let u_cart = FourVector {
        t: u_sph_next[0],
        x: vx,
        y: vy,
        z: vz,
    };
    assert_finite!(u_cart.t, "u_cart.t");
    assert_finite!(u_cart.x, "u_cart.x");
    assert_finite!(u_cart.y, "u_cart.y");
    assert_finite!(u_cart.z, "u_cart.z");

    (x_cart, u_cart)
}

// helper utility: cartesian_to_spherical kept for convenience
impl Metric {
    fn cartesian_to_spherical(x: &FourVector) -> (f64, f64, f64, f64) {
        let t = x.t;
        let r = (x.x * x.x + x.y * x.y + x.z * x.z).sqrt().max(1e-16);
        let theta = if r < 1e-12 {
            0.0
        } else {
            (x.z / r).clamp(-1.0, 1.0).acos()
        };
        let phi = if x.x.abs() < 1e-12 && x.y.abs() < 1e-12 {
            0.0
        } else {
            x.y.atan2(x.x)
        };
        (t, r, theta, phi)
    }
}

