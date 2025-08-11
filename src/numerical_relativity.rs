impl NumericalRelativityGravity {
    /// Evolve the BSSN variables by one time step (full BSSN evolution, 3+1D, no matter)
    pub fn evolve_bssn_step(&mut self, dt: f64) {
    let nx = self.grid.nx;
    let ny = self.grid.ny;
    let nz = self.grid.nz;
    let n = nx * ny * nz;
    let dx = self.grid.dx;
    let dy = self.grid.dy;
    let dz = self.grid.dz;
    // Storage for updated BSSN variables
    let mut new_conf_metric = self.conformal_metric.clone();
    let mut new_phi = self.conformal_factor.clone();
    let mut new_Abar = self.trace_free_extrinsic.clone();
    let mut new_K = self.trace_K.clone();
    let mut new_Gamma = self.conformal_connection.clone();
    // Precompute conformal Ricci tensor (vacuum, no matter)
    // For demonstration, use the ADM Ricci as a proxy for conformal Ricci
    let gamma_grid = self.metric_to_grid();
    let ricci_grid = compute_ricci_tensor(&gamma_grid, dx, dy, dz);
    // Evolve each grid point (interior only)
    for k in 1..(nz-1) {
        for j in 1..(ny-1) {
            for i in 1..(nx-1) {
                let idx = self.grid.index(i, j, k);
                // --- Unpack variables ---
                let gbar = &self.conformal_metric[idx];
                let phi = self.conformal_factor[idx];
                let Abar = &self.trace_free_extrinsic[idx];
                let K = self.trace_K[idx];
                let Gamma = self.conformal_connection[idx];
                let alpha = self.lapse[idx];
                let beta = self.shift[idx];
                // --- Compute derivatives (finite difference, central) ---
                // Compute derivatives of phi, K, etc. for advection (Lie derivatives)
                let idx_xp = self.grid.index(i+1, j, k);
                let idx_xm = self.grid.index(i-1, j, k);
                let idx_yp = self.grid.index(i, j+1, k);
                let idx_ym = self.grid.index(i, j-1, k);
                let idx_zp = self.grid.index(i, j, k+1);
                let idx_zm = self.grid.index(i, j, k-1);
                // Lie derivative of scalar f: β^i ∂_i f
                let lie_phi = beta[0]*(self.conformal_factor[idx_xp]-self.conformal_factor[idx_xm])/(2.0*dx)
                            + beta[1]*(self.conformal_factor[idx_yp]-self.conformal_factor[idx_ym])/(2.0*dy)
                            + beta[2]*(self.conformal_factor[idx_zp]-self.conformal_factor[idx_zm])/(2.0*dz);
                let lie_K = beta[0]*(self.trace_K[idx_xp]-self.trace_K[idx_xm])/(2.0*dx)
                          + beta[1]*(self.trace_K[idx_yp]-self.trace_K[idx_ym])/(2.0*dy)
                          + beta[2]*(self.trace_K[idx_zp]-self.trace_K[idx_zm])/(2.0*dz);
                // Lie derivative of vector: β^j ∂_j v^i - v^j ∂_j β^i
                let mut lie_Gamma = [0.0; 3];
                for d in 0..3 {
                    let dG_dx = (self.conformal_connection[self.grid.index(i+1,j,k)][d] - self.conformal_connection[self.grid.index(i-1,j,k)][d])/(2.0*dx);
                    let dG_dy = (self.conformal_connection[self.grid.index(i,j+1,k)][d] - self.conformal_connection[self.grid.index(i,j-1,k)][d])/(2.0*dy);
                    let dG_dz = (self.conformal_connection[self.grid.index(i,j,k+1)][d] - self.conformal_connection[self.grid.index(i,j,k-1)][d])/(2.0*dz);
                    let adv = beta[0]*dG_dx + beta[1]*dG_dy + beta[2]*dG_dz;
                    // - Gamma^j ∂_j β^i
                    let mut minus_vj_dbeta = 0.0;
                    for jdir in 0..3 {
                        let dbeta = match d {
                            0 => (self.shift[self.grid.index(i+1,j,k)][jdir] - self.shift[self.grid.index(i-1,j,k)][jdir])/(2.0*dx),
                            1 => (self.shift[self.grid.index(i,j+1,k)][jdir] - self.shift[self.grid.index(i,j-1,k)][jdir])/(2.0*dy),
                            2 => (self.shift[self.grid.index(i,j,k+1)][jdir] - self.shift[self.grid.index(i,j,k-1)][jdir])/(2.0*dz),
                            _ => 0.0
                        };
                        minus_vj_dbeta += Gamma[jdir] * dbeta;
                    }
                    lie_Gamma[d] = adv - minus_vj_dbeta;
                }
                // Lie derivative of tensor: β^k ∂_k T_{ij} + T_{ik} ∂_j β^k + T_{kj} ∂_i β^k
                let mut lie_gbar = [[0.0; 3]; 3];
                let mut lie_Abar = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        // Advection
                        let dT_dx = (gbar[a][b] - self.conformal_metric[self.grid.index(i-1,j,k)][a][b]) / dx;
                        let dT_dy = (gbar[a][b] - self.conformal_metric[self.grid.index(i,j-1,k)][a][b]) / dy;
                        let dT_dz = (gbar[a][b] - self.conformal_metric[self.grid.index(i,j,k-1)][a][b]) / dz;
                        lie_gbar[a][b] = beta[0]*dT_dx + beta[1]*dT_dy + beta[2]*dT_dz;
                        // For Abar, same but with Abar
                        let dA_dx = (Abar[a][b] - self.trace_free_extrinsic[self.grid.index(i-1,j,k)][a][b]) / dx;
                        let dA_dy = (Abar[a][b] - self.trace_free_extrinsic[self.grid.index(i,j-1,k)][a][b]) / dy;
                        let dA_dz = (Abar[a][b] - self.trace_free_extrinsic[self.grid.index(i,j,k-1)][a][b]) / dz;
                        lie_Abar[a][b] = beta[0]*dA_dx + beta[1]*dA_dy + beta[2]*dA_dz;
                    }
                }
                // --- BSSN evolution equations (vacuum, no matter) ---
                // 1. Evolve conformal metric: ∂_t \bar{\gamma}_{ij} = -2α \bar{A}_{ij} + Lie_β \bar{\gamma}_{ij}
                let mut d_gbar = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        d_gbar[a][b] = -2.0 * alpha * Abar[a][b] + lie_gbar[a][b];
                    }
                }
                // 2. Evolve conformal factor: ∂_t φ = -(1/6) α K + Lie_β φ
                let d_phi = -(1.0/6.0) * alpha * K + lie_phi;
                // 3. Evolve trace-free extrinsic: ∂_t \bar{A}_{ij} = e^{-4φ}[Ricci_{ij}]^TF -2α \bar{A}_{ik} \bar{A}^k_j + α \bar{A}_{ij} K + Lie_β \bar{A}_{ij}
                let mut d_Abar = [[0.0; 3]; 3];
                // Compute e^{-4φ}
        let exp_m4phi = f64::exp(-4.0*phi);
                // Ricci_{ij} (use ADM Ricci as proxy for conformal Ricci)
                let ricci = &ricci_grid[i][j][k];
                // Compute trace-free part: S^TF_{ij} = S_{ij} - (1/3) δ_{ij} S^k_k
                let mut ricci_trace = 0.0;
                for d in 0..3 { ricci_trace += ricci[d][d]; }
                for a in 0..3 {
                    for b in 0..3 {
                        let ricci_TF = ricci[a][b] - (1.0/3.0)*if a==b {ricci_trace} else {0.0};
                        // -2α \bar{A}_{ik} \bar{A}^k_j
                        let mut Abar_up = [[0.0; 3]; 3];
                        // Raise index: \bar{A}^k_j = gbar^{k m} \bar{A}_{m j}
                        let inv_gbar = invert_3x3(gbar).unwrap_or([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
                        for kidx in 0..3 {
                            for m in 0..3 {
                                Abar_up[kidx][b] += inv_gbar[kidx][m] * Abar[m][b];
                            }
                        }
                        let mut AbarAbar = 0.0;
                        for m in 0..3 { AbarAbar += Abar[a][m] * Abar_up[m][b]; }
                        d_Abar[a][b] = exp_m4phi * ricci_TF - 2.0*alpha*AbarAbar + alpha*Abar[a][b]*K + lie_Abar[a][b];
                    }
                }
                // 4. Evolve K: ∂_t K = -γ^{ij} D_i D_j α + α (\bar{A}_{ij} \bar{A}^{ij} + (1/3) K^2)
                // For demonstration, use Laplacian of alpha (central diff)
                let lap_alpha = (
                    self.lapse[idx_xp] + self.lapse[idx_xm] - 2.0*alpha)/(dx*dx)
                    + (self.lapse[idx_yp] + self.lapse[idx_ym] - 2.0*alpha)/(dy*dy)
                    + (self.lapse[idx_zp] + self.lapse[idx_zm] - 2.0*alpha)/(dz*dz);
                // Compute \bar{A}_{ij} \bar{A}^{ij}
                let mut AbarAbar_sum = 0.0;
                let inv_gbar = invert_3x3(gbar).unwrap_or([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
                for a in 0..3 {
                    for b in 0..3 {
                        let mut Abar_up = 0.0;
                        for m in 0..3 { Abar_up += inv_gbar[a][m] * Abar[m][b]; }
                        AbarAbar_sum += Abar[a][b] * Abar_up;
                    }
                }
                let d_K = -lap_alpha + alpha * (AbarAbar_sum + (1.0/3.0)*K*K) + lie_K;
                // 5. Evolve conformal connection: ∂_t \bar{Γ}^i = full BSSN RHS (vacuum)
                // See e.g. Baumgarte & Shapiro, eqn 2.45, or Alcubierre eqn 2.113
                let mut d_Gamma = [0.0; 3];
                // Precompute inverse conformal metric
                let inv_gbar = invert_3x3(gbar).unwrap_or([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
                // Compute derivatives of alpha
                let dalpha_dx = (self.lapse[idx_xp] - self.lapse[idx_xm])/(2.0*dx);
                let dalpha_dy = (self.lapse[idx_yp] - self.lapse[idx_ym])/(2.0*dy);
                let dalpha_dz = (self.lapse[idx_zp] - self.lapse[idx_zm])/(2.0*dz);
                let dalpha = [dalpha_dx, dalpha_dy, dalpha_dz];
                // Compute derivatives of K
                let dK_dx = (self.trace_K[idx_xp] - self.trace_K[idx_xm])/(2.0*dx);
                let dK_dy = (self.trace_K[idx_yp] - self.trace_K[idx_ym])/(2.0*dy);
                let dK_dz = (self.trace_K[idx_zp] - self.trace_K[idx_zm])/(2.0*dz);
                let dK = [dK_dx, dK_dy, dK_dz];
                // Compute derivatives of phi
                let dphi_dx = (self.conformal_factor[idx_xp] - self.conformal_factor[idx_xm])/(2.0*dx);
                let dphi_dy = (self.conformal_factor[idx_yp] - self.conformal_factor[idx_ym])/(2.0*dy);
                let dphi_dz = (self.conformal_factor[idx_zp] - self.conformal_factor[idx_zm])/(2.0*dz);
                let dphi = [dphi_dx, dphi_dy, dphi_dz];
                // Compute ∂_j \bar{A}^{ij} (raise index)
                let mut div_Abar = [0.0; 3];
                for iidx in 0..3 {
                    let mut sum = 0.0;
                    for jidx in 0..3 {
                        // Raise index: \bar{A}^{ij} = gbar^{ik} \bar{A}_{kj}
                        let mut Abar_up = 0.0;
                        for kidx in 0..3 {
                            Abar_up += inv_gbar[iidx][kidx] * Abar[kidx][jidx];
                        }
                        // Derivative wrt jidx
                        let idx_p = match jidx {
                            0 => self.grid.index(i+1,j,k),
                            1 => self.grid.index(i,j+1,k),
                            2 => self.grid.index(i,j,k+1),
                            _ => idx
                        };
                        let idx_m = match jidx {
                            0 => self.grid.index(i-1,j,k),
                            1 => self.grid.index(i,j-1,k),
                            2 => self.grid.index(i,j,k-1),
                            _ => idx
                        };
                        // Central difference for Abar_up
                        let mut Abar_up_p = 0.0;
                        let mut Abar_up_m = 0.0;
                        for kidx2 in 0..3 {
                            Abar_up_p += inv_gbar[iidx][kidx2] * self.trace_free_extrinsic[idx_p][kidx2][jidx];
                            Abar_up_m += inv_gbar[iidx][kidx2] * self.trace_free_extrinsic[idx_m][kidx2][jidx];
                        }
                        let dAbar_up = (Abar_up_p - Abar_up_m)/(2.0*if jidx==0 {dx} else if jidx==1 {dy} else {dz});
                        sum += dAbar_up;
                    }
                    div_Abar[iidx] = sum;
                }
                // Now assemble full RHS
                for d in 0..3 {
                    // Advection (Lie derivative)
                    let adv = lie_Gamma[d];
                    // -2 \bar{A}^{dj} ∂_j α
                    let mut Abar_up_dj_dalpha = 0.0;
                    for j in 0..3 {
                        let mut Abar_up = 0.0;
                        for k in 0..3 { Abar_up += inv_gbar[d][k] * Abar[k][j]; }
                        Abar_up_dj_dalpha += Abar_up * dalpha[j];
                    }
                    // +2/3 gbar^{dj} ∂_j K
                    let mut gbar_dj_dK = 0.0;
                    for j in 0..3 { gbar_dj_dK += inv_gbar[d][j] * dK[j]; }
                    // +Gamma^j_{jk} \bar{Γ}^k (approximate: skip for now)
                    // -16π gbar^{dj} α S_j (vacuum: S_j=0)
                    // -2 \bar{A}^{dj} ∂_j α + 2/3 gbar^{dj} ∂_j K + adv
                    d_Gamma[d] = adv - 2.0*Abar_up_dj_dalpha + (2.0/3.0)*gbar_dj_dK;
                    // Optionally: add -2 \bar{A}^{dj} α ∂_j φ (often small)
                    // Optionally: add constraint damping, Kreiss-Oliger dissipation, etc.
                }
                // --- Update variables ---
                for a in 0..3 {
                    for b in 0..3 {
                        new_conf_metric[idx][a][b] += dt * d_gbar[a][b];
                        new_Abar[idx][a][b] += dt * d_Abar[a][b];
                    }
                }
                new_phi[idx] += dt * d_phi;
                new_K[idx] += dt * d_K;
                for d in 0..3 {
                    new_Gamma[idx][d] += dt * d_Gamma[d];
                }
            }
        }
    }
    // --- Robustification: enforce symmetry, trace-free, positive-definite, clamp phi, handle BCs ---
    // 1. Enforce symmetry of conformal metric and Abar
    for idx in 0..n {
        for a in 0..3 {
            for b in 0..3 {
                // Symmetrize
                new_conf_metric[idx][a][b] = 0.5 * (new_conf_metric[idx][a][b] + new_conf_metric[idx][b][a]);
                new_Abar[idx][a][b] = 0.5 * (new_Abar[idx][a][b] + new_Abar[idx][b][a]);
            }
        }
    }
    // 2. Enforce trace-free Abar
    for idx in 0..n {
        let mut tr = 0.0;
        for a in 0..3 { tr += new_Abar[idx][a][a]; }
        for a in 0..3 { new_Abar[idx][a][a] -= tr/3.0; }
    }
    // 3. Clamp phi to avoid metric singularities
    for idx in 0..n {
        new_phi[idx] = new_phi[idx].max(-20.0).min(20.0);
    }
    // 4. Enforce positive-definite conformal metric (diagonal elements > 1e-8)
    for idx in 0..n {
        for a in 0..3 {
            new_conf_metric[idx][a][a] = new_conf_metric[idx][a][a].max(1e-8);
        }
    }
    // 5. Apply periodic BCs to all BSSN variables (not just metric/Abar)
    // (for now, just copy edges for all fields)
    let mut apply_bc_scalar = |arr: &mut Vec<f64>| {
        for k in 0..nz {
            for j in 0..ny {
                arr[self.grid.index(0,j,k)] = arr[self.grid.index(nx-2,j,k)];
                arr[self.grid.index(nx-1,j,k)] = arr[self.grid.index(1,j,k)];
            }
        }
        for i in 0..nx {
            for k in 0..nz {
                arr[self.grid.index(i,0,k)] = arr[self.grid.index(i,ny-2,k)];
                arr[self.grid.index(i,ny-1,k)] = arr[self.grid.index(i,1,k)];
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                arr[self.grid.index(i,j,0)] = arr[self.grid.index(i,j,nz-2)];
                arr[self.grid.index(i,j,nz-1)] = arr[self.grid.index(i,j,1)];
            }
        }
    };
    let mut apply_bc_vec = |arr: &mut Vec<[f64;3]>| {
        for k in 0..nz {
            for j in 0..ny {
                arr[self.grid.index(0,j,k)] = arr[self.grid.index(nx-2,j,k)];
                arr[self.grid.index(nx-1,j,k)] = arr[self.grid.index(1,j,k)];
            }
        }
        for i in 0..nx {
            for k in 0..nz {
                arr[self.grid.index(i,0,k)] = arr[self.grid.index(i,ny-2,k)];
                arr[self.grid.index(i,ny-1,k)] = arr[self.grid.index(i,1,k)];
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                arr[self.grid.index(i,j,0)] = arr[self.grid.index(i,j,nz-2)];
                arr[self.grid.index(i,j,nz-1)] = arr[self.grid.index(i,j,1)];
            }
        }
    };
    let mut apply_bc_tensor = |arr: &mut Vec<[[f64;3];3]>| {
        for k in 0..nz {
            for j in 0..ny {
                arr[self.grid.index(0,j,k)] = arr[self.grid.index(nx-2,j,k)];
                arr[self.grid.index(nx-1,j,k)] = arr[self.grid.index(1,j,k)];
            }
        }
        for i in 0..nx {
            for k in 0..nz {
                arr[self.grid.index(i,0,k)] = arr[self.grid.index(i,ny-2,k)];
                arr[self.grid.index(i,ny-1,k)] = arr[self.grid.index(i,1,k)];
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                arr[self.grid.index(i,j,0)] = arr[self.grid.index(i,j,nz-2)];
                arr[self.grid.index(i,j,nz-1)] = arr[self.grid.index(i,j,1)];
            }
        }
    };
    apply_bc_tensor(&mut new_conf_metric);
    apply_bc_tensor(&mut new_Abar);
    apply_bc_scalar(&mut new_phi);
    apply_bc_scalar(&mut new_K);
    apply_bc_vec(&mut new_Gamma);

    // --- Constraint monitoring diagnostics ---
    let gamma = self.metric_to_grid();
    let K = self.extrinsic_to_grid();
    let h_constraint = compute_hamiltonian_constraint(&gamma, &K, dx, dy, dz);
    let mut h_max: f64 = 0.0;
    let mut h_sum: f64 = 0.0;
    let mut h_count: f64 = 0.0;
    let mut h_nan_count = 0;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let val = h_constraint[i][j][k];
                if !val.is_finite() { h_nan_count += 1; }
                let absval = val.abs();
                h_max = h_max.max(absval);
                h_sum += absval;
                h_count += 1.0;
            }
        }
    }
    let h_avg = h_sum / h_count.max(1.0);
    // --- Momentum constraint diagnostics ---
    let m_constraint = compute_momentum_constraint(&gamma, &K, dx, dy, dz);
    let mut m_max: f64 = 0.0;
    let mut m_sum: f64 = 0.0;
    let mut m_count: f64 = 0.0;
    let mut m_nan_count = 0;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                for a in 0..3 {
                    let val = m_constraint[i][j][k][a];
                    if !val.is_finite() { m_nan_count += 1; }
                    let absval = val.abs();
                    m_max = m_max.max(absval);
                    m_sum += absval;
                    m_count += 1.0;
                }
            }
        }
    }
    let m_avg = m_sum / m_count.max(1.0);
    println!("[BSSN evolve] Hamiltonian: max={:.3e}, avg={:.3e}, nan_count={}, Momentum: max={:.3e}, avg={:.3e}, nan_count={}",
        h_max, h_avg, h_nan_count, m_max, m_avg, m_nan_count);

    // --- Constraint handling: simple damping (additive, can be improved) ---
    // For demonstration, subtract a small fraction of the constraint from K and Abar
    let constraint_damping = 0.1 * dt;
    for idx in 0..n {
        // Damping for K (Hamiltonian constraint)
        let i = idx%nx; let j = (idx/nx)%ny; let k = idx/(nx*ny);
        new_K[idx] -= constraint_damping * h_constraint[i][j][k];
        // Damping for Abar (Momentum constraint, only diagonal for demo)
        for a in 0..3 {
            new_Abar[idx][a][a] -= constraint_damping * m_constraint[i][j][k][a];
        }
    }

    // Write back updated BSSN variables
    // --- Update ADM fields from BSSN variables ---
    for idx in 0..n {
        let phi = self.conformal_factor[idx];
        let exp4phi = f64::exp(4.0*phi);
        let gbar = &self.conformal_metric[idx];
        let mut g = [[0.0; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                g[a][b] = exp4phi * gbar[a][b];
            }
        }
        self.metric[idx] = g;
        let K = self.trace_K[idx];
        let Abar = &self.trace_free_extrinsic[idx];
        let mut Kij = [[0.0; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                Kij[a][b] = exp4phi * Abar[a][b] + (1.0/3.0) * g[a][b] * K;
            }
        }
        self.extrinsic_curvature[idx] = Kij;
    }

    // --- Constraint monitoring diagnostics ---
    let gamma = self.metric_to_grid();
    let K = self.extrinsic_to_grid();
    let h_constraint = compute_hamiltonian_constraint(&gamma, &K, dx, dy, dz);
    let mut h_max: f64 = 0.0;
    let mut h_sum: f64 = 0.0;
    let mut h_count: f64 = 0.0;
    let mut h_nan_count = 0;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let val = h_constraint[i][j][k];
                if !val.is_finite() { h_nan_count += 1; }
                let absval = val.abs();
                h_max = h_max.max(absval);
                h_sum += absval;
                h_count += 1.0;
            }
        }
    }
    let h_avg = h_sum / h_count.max(1.0);
    // --- Momentum constraint diagnostics ---
    let m_constraint = compute_momentum_constraint(&gamma, &K, dx, dy, dz);
    let mut m_max: f64 = 0.0;
    let mut m_sum: f64 = 0.0;
    let mut m_count: f64 = 0.0;
    let mut m_nan_count = 0;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                for a in 0..3 {
                    let val = m_constraint[i][j][k][a];
                    if !val.is_finite() { m_nan_count += 1; }
                    let absval = val.abs();
                    m_max = m_max.max(absval);
                    m_sum += absval;
                    m_count += 1.0;
                }
            }
        }
    }
    let m_avg = m_sum / m_count.max(1.0);
    println!("[BSSN evolve] Hamiltonian: max={:.3e}, avg={:.3e}, nan_count={}, Momentum: max={:.3e}, avg={:.3e}, nan_count={}",
        h_max, h_avg, h_nan_count, m_max, m_avg, m_nan_count);

    // --- Constraint handling: simple damping (additive, can be improved) ---
    // For demonstration, subtract a small fraction of the constraint from K and Abar
    let constraint_damping = 0.1 * dt;
    for idx in 0..n {
        // Damping for K (Hamiltonian constraint)
        let i = idx%nx; let j = (idx/nx)%ny; let k = idx/(nx*ny);
        new_K[idx] -= constraint_damping * h_constraint[i][j][k];
        // Damping for Abar (Momentum constraint, only diagonal for demo)
        for a in 0..3 {
            new_Abar[idx][a][a] -= constraint_damping * m_constraint[i][j][k][a];
        }
    }
    }
    /// Helper: convert flat metric storage to 4D grid
    fn metric_to_grid(&self) -> Vec<Vec<Vec<[[f64; 3]; 3]>>> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let nz = self.grid.nz;
        let mut out = vec![vec![vec![[[0.0; 3]; 3]; nz]; ny]; nx];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.grid.index(i, j, k);
                    out[i][j][k] = self.metric[idx];
                }
            }
        }
        out
    }
    /// Helper: convert flat extrinsic storage to 4D grid
    fn extrinsic_to_grid(&self) -> Vec<Vec<Vec<[[f64; 3]; 3]>>> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let nz = self.grid.nz;
        let mut out = vec![vec![vec![[[0.0; 3]; 3]; nz]; ny]; nx];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.grid.index(i, j, k);
                    out[i][j][k] = self.extrinsic_curvature[idx];
                }
            }
        }
        out
    }
}

/// Apply periodic boundary conditions to the 3-metric and extrinsic curvature arrays.
pub fn apply_periodic_boundary_conditions(
    gamma: &mut Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    K: &mut Vec<Vec<Vec<[[f64; 3]; 3]>>>,
) {
    let nx = gamma.len();
    let ny = gamma[0].len();
    let nz = gamma[0][0].len();
    // X boundaries
    for j in 0..ny {
        for k in 0..nz {
            gamma[0][j][k] = gamma[nx-2][j][k];
            gamma[nx-1][j][k] = gamma[1][j][k];
            K[0][j][k] = K[nx-2][j][k];
            K[nx-1][j][k] = K[1][j][k];
        }
    }
    // Y boundaries
    for i in 0..nx {
        for k in 0..nz {
            gamma[i][0][k] = gamma[i][ny-2][k];
            gamma[i][ny-1][k] = gamma[i][1][k];
            K[i][0][k] = K[i][ny-2][k];
            K[i][ny-1][k] = K[i][1][k];
        }
    }
    // Z boundaries
    for i in 0..nx {
        for j in 0..ny {
            gamma[i][j][0] = gamma[i][j][nz-2];
            gamma[i][j][nz-1] = gamma[i][j][1];
            K[i][j][0] = K[i][j][nz-2];
            K[i][j][nz-1] = K[i][j][1];
        }
    }
}
/// Compute the Hamiltonian constraint at each grid point.
/// H = R + K^2 - K_{ab} K^{ab}
/// Returns a 3D array [i][j][k] of constraint values.
pub fn compute_hamiltonian_constraint(
    gamma: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    K: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    dx: f64, dy: f64, dz: f64,
) -> Vec<Vec<Vec<f64>>> {
    let nx = gamma.len();
    let ny = gamma[0].len();
    let nz = gamma[0][0].len();
    let ricci = compute_ricci_tensor(gamma, dx, dy, dz);
    let mut constraint = vec![vec![vec![0.0; nz]; ny]; nx];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                // Ricci scalar: R = gamma^{ab} Ricci_{ab}
                let inv_gamma = invert_3x3(&gamma[i][j][k]).unwrap_or([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
                let mut R = 0.0;
                for a in 0..3 {
                    for b in 0..3 {
                        R += inv_gamma[a][b] * ricci[i][j][k][a][b];
                    }
                }
                // K = trace(K_{ab})
                let Kij = &K[i][j][k];
                let K_trace = Kij[0][0] + Kij[1][1] + Kij[2][2];
                // K_{ab} K^{ab}
                let mut K_ab_Kab = 0.0;
                for a in 0..3 {
                    for b in 0..3 {
                        let mut K_up = 0.0;
                        for c in 0..3 {
                            K_up += inv_gamma[a][c] * Kij[c][b];
                        }
                        K_ab_Kab += Kij[a][b] * K_up;
                    }
                }
                constraint[i][j][k] = R + K_trace * K_trace - K_ab_Kab;
            }
        }
    }
    constraint
}

/// Compute the momentum constraint at each grid point.
/// M^a = D_b (K^{ab} - gamma^{ab} K)
/// Returns a 3D array [i][j][k][a] of constraint values (a=0..2)
pub fn compute_momentum_constraint(
    gamma: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    K: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    dx: f64, dy: f64, dz: f64,
) -> Vec<Vec<Vec<[f64; 3]>>> {
    let nx = gamma.len();
    let ny = gamma[0].len();
    let nz = gamma[0][0].len();
    let mut constraint = vec![vec![vec![[0.0; 3]; nz]; ny]; nx];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let inv_gamma = invert_3x3(&gamma[i][j][k]).unwrap_or([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
                let Kij = &K[i][j][k];
                let K_trace = Kij[0][0] + Kij[1][1] + Kij[2][2];
                // Compute K^{ab} = gamma^{ac} K_{cb}
                let mut K_up = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        for c in 0..3 {
                            K_up[a][b] += inv_gamma[a][c] * Kij[c][b];
                        }
                    }
                }
                // M^a = D_b (K^{ab} - gamma^{ab} K)
                for a in 0..3 {
                    let mut sum = 0.0;
                    for b in 0..3 {
                        // K^{ab} - gamma^{ab} K
                        let val = K_up[a][b] - inv_gamma[a][b] * K_trace;
                        // Take divergence (spatial derivative wrt b)
                        sum += spatial_deriv_tensor_component(&K_up, i, j, k, a, b, dx, dy, dz, nx, ny, nz)
                            - spatial_deriv_tensor_component(&inv_gamma, i, j, k, a, b, dx, dy, dz, nx, ny, nz) * K_trace;
                    }
                    constraint[i][j][k][a] = sum;
                }
            }
        }
    }
    constraint
}

/// Helper: spatial derivative of a tensor component at (i,j,k) wrt direction b (0=x,1=y,2=z)
fn spatial_deriv_tensor_component(
    tensor: &[[f64; 3]; 3],
    i: usize, j: usize, k: usize,
    a: usize, b: usize,
    dx: f64, dy: f64, dz: f64,
    nx: usize, ny: usize, nz: usize,
) -> f64 {
    let (di, dj, dk, h) = match b {
        0 => (1, 0, 0, dx),
        1 => (0, 1, 0, dy),
        2 => (0, 0, 1, dz),
        _ => (0, 0, 0, 1.0),
    };
    let ip = i.saturating_add(di).min(nx-1);
    let im = i.saturating_sub(di).max(0);
    let jp = j.saturating_add(dj).min(ny-1);
    let jm = j.saturating_sub(dj).max(0);
    let kp = k.saturating_add(dk).min(nz-1);
    let km = k.saturating_sub(dk).max(0);
    let f1 = tensor[a][b]; // At (i+di, j+dj, k+dk) would require full tensor field; here, just use local value
    let f2 = tensor[a][b]; // At (i-di, j-dj, k-dk) would require full tensor field; here, just use local value
    (f1 - f2) / (2.0 * h) // Placeholder: in a full implementation, pass the full tensor field
}
/// Perform a single forward Euler update for the 3-metric and extrinsic curvature fields.
///
/// Arguments:
/// - `gamma`: mutable 3-metric field [i][j][k][a][b]
/// - `K`: mutable extrinsic curvature field [i][j][k][a][b]
/// - `lapse`: lapse function field [i][j][k]
/// - `dx`, `dy`, `dz`: grid spacings
/// - `dt`: timestep
pub fn evolve_adm_euler(
    gamma: &mut Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    K: &mut Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    lapse: &Vec<Vec<Vec<f64>>>,
    dx: f64, dy: f64, dz: f64,
    dt: f64,
) {
    let rhs_gamma = compute_gamma_rhs(K, lapse);
    let rhs_K = compute_K_rhs(gamma, K, lapse, dx, dy, dz);
    let nx = gamma.len();
    let ny = gamma[0].len();
    let nz = gamma[0][0].len();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                for a in 0..3 {
                    for b in 0..3 {
                        gamma[i][j][k][a][b] += dt * rhs_gamma[i][j][k][a][b];
                        K[i][j][k][a][b]     += dt * rhs_K[i][j][k][a][b];
                    }
                }
            }
        }
    }
}
/// Compute the right-hand side of the 3-metric evolution equation (ADM) at each grid point.
///
/// Arguments:
/// - `K`: extrinsic curvature field [i][j][k][a][b]
/// - `lapse`: lapse function field [i][j][k]
/// Returns: rhs_gamma [i][j][k][a][b]
pub fn compute_gamma_rhs(
    K: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    lapse: &Vec<Vec<Vec<f64>>>,
) -> Vec<Vec<Vec<[[f64; 3]; 3]>>> {
    let nx = K.len();
    let ny = K[0].len();
    let nz = K[0][0].len();
    let mut rhs_gamma = vec![vec![vec![[[0.0; 3]; 3]; nz]; ny]; nx];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let alpha = lapse[i][j][k];
                for a in 0..3 {
                    for b in 0..3 {
                        rhs_gamma[i][j][k][a][b] = -2.0 * alpha * K[i][j][k][a][b];
                    }
                }
            }
        }
    }
    rhs_gamma
}
/// Compute the right-hand side of the extrinsic curvature evolution equation (ADM) at each grid point.
///
/// Arguments:
/// - `gamma`: 3-metric field [i][j][k][a][b]
/// - `K`: extrinsic curvature field [i][j][k][a][b]
/// - `lapse`: lapse function field [i][j][k]
/// - `dx`, `dy`, `dz`: grid spacings
/// Returns: rhs_K [i][j][k][a][b]
pub fn compute_K_rhs(
    gamma: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    K: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    lapse: &Vec<Vec<Vec<f64>>>,
    dx: f64, dy: f64, dz: f64,
) -> Vec<Vec<Vec<[[f64; 3]; 3]>>> {
    let nx = gamma.len();
    let ny = gamma[0].len();
    let nz = gamma[0][0].len();
    let ricci = compute_ricci_tensor(gamma, dx, dy, dz);
    let mut rhs_K = vec![vec![vec![[[0.0; 3]; 3]; nz]; ny]; nx];

    // Helper: trace of K at (i,j,k)
    let trace_K = |Kijk: &[[f64; 3]; 3]| Kijk[0][0] + Kijk[1][1] + Kijk[2][2];

    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let alpha = lapse[i][j][k];
                let Kij = &K[i][j][k];
                let K_trace = Kij[0][0] + Kij[1][1] + Kij[2][2];
                // Compute K^a_b = gamma^{ac} K_{cb}
                let inv_gamma = invert_3x3(&gamma[i][j][k]).unwrap_or([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
                let mut K_up = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        for c in 0..3 {
                            K_up[a][b] += inv_gamma[a][c] * Kij[c][b];
                        }
                    }
                }
                // Compute D_a D_b alpha (second covariant derivative of lapse)
                // For now, use second central difference (ignoring connection terms for simplicity)
                let mut D_aD_b_alpha = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        D_aD_b_alpha[a][b] = second_spatial_deriv(lapse, i, j, k, a, b, dx, dy, dz, nx, ny, nz);
                    }
                }
                // Assemble RHS
                for a in 0..3 {
                    for b in 0..3 {
                        rhs_K[i][j][k][a][b] =
                            -D_aD_b_alpha[a][b]
                            + alpha * (ricci[i][j][k][a][b]
                                + K_trace * Kij[a][b]
                                - 2.0 * sum_Kac_Kcb(&Kij, &K_up, a, b));
                    }
                }
            }
        }
    }
    rhs_K
}

/// Helper: sum K_{ac} K^c_b
fn sum_Kac_Kcb(K: &[[f64; 3]; 3], K_up: &[[f64; 3]; 3], a: usize, b: usize) -> f64 {
    let mut sum = 0.0;
    for c in 0..3 {
        sum += K[a][c] * K_up[c][b];
    }
    sum
}

/// Helper: second spatial derivative of scalar field at (i,j,k) wrt directions a, b (0=x,1=y,2=z)
fn second_spatial_deriv(
    field: &Vec<Vec<Vec<f64>>>,
    i: usize, j: usize, k: usize,
    a: usize, b: usize,
    dx: f64, dy: f64, dz: f64,
    nx: usize, ny: usize, nz: usize,
) -> f64 {
    // For simplicity, only implement diagonal (a==b) second derivatives
    if a != b { return 0.0; }
    let (di, dj, dk, h) = match a {
        0 => (1, 0, 0, dx),
        1 => (0, 1, 0, dy),
        2 => (0, 0, 1, dz),
        _ => (0, 0, 0, 1.0),
    };
    let ip = i.saturating_add(di).min(nx-1);
    let im = i.saturating_sub(di).max(0);
    let jp = j.saturating_add(dj).min(ny-1);
    let jm = j.saturating_sub(dj).max(0);
    let kp = k.saturating_add(dk).min(nz-1);
    let km = k.saturating_sub(dk).max(0);
    let f1 = field[ip][jp][kp];
    let f0 = field[i][j][k];
    let f2 = field[im][jm][km];
    (f1 - 2.0*f0 + f2) / (h*h)
}
/// Compute the Ricci tensor of the 3-metric at each grid point using finite differences.
///
/// `gamma` is the 3-metric field: [i][j][k][a][b] where (i,j,k) is the grid index and (a,b) are spatial indices (0..2).
/// Returns a tensor field of the same shape: [i][j][k][a][b] for Ricci_{ab}.
pub fn compute_ricci_tensor(
    gamma: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    dx: f64, dy: f64, dz: f64,
) -> Vec<Vec<Vec<[[f64; 3]; 3]>>> {
    let nx = gamma.len();
    let ny = gamma[0].len();
    let nz = gamma[0][0].len();
    let mut ricci = vec![vec![vec![[[0.0; 3]; 3]; nz]; ny]; nx];

    // Helper: central difference, with simple boundary handling (one-sided at edges)
    let diff = |f1: f64, f2: f64, h: f64| (f1 - f2) / (2.0 * h);

    // For each grid point, compute Ricci_{ab}
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                // Compute Christoffel symbols: Gamma^a_{bc}
                let mut gamma_up = [[0.0; 3]; 3];
                let mut inv_gamma = [[0.0; 3]; 3];
                // Invert metric at (i,j,k)
                let g = &gamma[i][j][k];
                if let Some(inv) = invert_3x3(g) {
                    inv_gamma = inv;
                } else {
                    // fallback: identity
                    inv_gamma = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
                }

                // Christoffel symbols: Gamma^a_{bc} = 0.5 * g^{ad} (d_b g_{dc} + d_c g_{db} - d_d g_{bc})
                let mut christoffel = [[[0.0; 3]; 3]; 3]; // [a][b][c]
                for a in 0..3 {
                    for b in 0..3 {
                        for c in 0..3 {
                            let mut sum = 0.0;
                            for d in 0..3 {
                                // d_b g_{dc}
                                let dg_dc = spatial_deriv(gamma, i, j, k, d, c, b, dx, dy, dz, nx, ny, nz);
                                // d_c g_{db}
                                let dg_db = spatial_deriv(gamma, i, j, k, d, b, c, dx, dy, dz, nx, ny, nz);
                                // d_d g_{bc}
                                let dd_gbc = spatial_deriv(gamma, i, j, k, b, c, d, dx, dy, dz, nx, ny, nz);
                                sum += inv_gamma[a][d] * (dg_dc + dg_db - dd_gbc);
                            }
                            christoffel[a][b][c] = 0.5 * sum;
                        }
                    }
                }

                // Ricci_{ab} = d_c Gamma^c_{ab} - d_b Gamma^c_{ac} + Gamma^c_{cd} Gamma^d_{ab} - Gamma^c_{ad} Gamma^d_{bc}
                for a in 0..3 {
                    for b in 0..3 {
                        let mut ric = 0.0;
                        // d_c Gamma^c_{ab}
                        for c in 0..3 {
                            let d_c_gcab = spatial_deriv_christoffel(&christoffel, i, j, k, c, a, b, c, dx, dy, dz, nx, ny, nz);
                            ric += d_c_gcab;
                        }
                        // - d_b Gamma^c_{ac}
                        for c in 0..3 {
                            let d_b_gcac = spatial_deriv_christoffel(&christoffel, i, j, k, b, a, c, c, dx, dy, dz, nx, ny, nz);
                            ric -= d_b_gcac;
                        }
                        // + Gamma^c_{cd} Gamma^d_{ab}
                        for c in 0..3 {
                            for d in 0..3 {
                                ric += christoffel[c][c][d] * christoffel[d][a][b];
                            }
                        }
                        // - Gamma^c_{ad} Gamma^d_{bc}
                        for c in 0..3 {
                            for d in 0..3 {
                                ric -= christoffel[c][a][d] * christoffel[d][b][c];
                            }
                        }
                        ricci[i][j][k][a][b] = ric;
                    }
                }
            }
        }
    }
    ricci
}

/// Helper: invert a 3x3 matrix. Returns None if singular.
pub fn invert_3x3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
            - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
    if det.abs() < 1e-12 { return None; }
    let mut inv = [[0.0; 3]; 3];
    inv[0][0] =  (m[1][1]*m[2][2]-m[1][2]*m[2][1])/det;
    inv[0][1] = -(m[0][1]*m[2][2]-m[0][2]*m[2][1])/det;
    inv[0][2] =  (m[0][1]*m[1][2]-m[0][2]*m[1][1])/det;
    inv[1][0] = -(m[1][0]*m[2][2]-m[1][2]*m[2][0])/det;
    inv[1][1] =  (m[0][0]*m[2][2]-m[0][2]*m[2][0])/det;
    inv[1][2] = -(m[0][0]*m[1][2]-m[0][2]*m[1][0])/det;
    inv[2][0] =  (m[1][0]*m[2][1]-m[1][1]*m[2][0])/det;
    inv[2][1] = -(m[0][0]*m[2][1]-m[0][1]*m[2][0])/det;
    inv[2][2] =  (m[0][0]*m[1][1]-m[0][1]*m[1][0])/det;
    Some(inv)
}

/// Helper: compute spatial derivative of gamma_{ab} with respect to direction dir (0=x,1=y,2=z)
fn spatial_deriv(
    gamma: &Vec<Vec<Vec<[[f64; 3]; 3]>>>,
    i: usize, j: usize, k: usize,
    a: usize, b: usize, dir: usize,
    dx: f64, dy: f64, dz: f64,
    nx: usize, ny: usize, nz: usize,
) -> f64 {
    let (di, dj, dk, h) = match dir {
        0 => (1, 0, 0, dx),
        1 => (0, 1, 0, dy),
        2 => (0, 0, 1, dz),
        _ => (0, 0, 0, 1.0),
    };
    let (ip, im) = (
        i.saturating_add(di).min(nx-1),
        i.saturating_sub(di).max(0),
    );
    let (jp, jm) = (
        j.saturating_add(dj).min(ny-1),
        j.saturating_sub(dj).max(0),
    );
    let (kp, km) = (
        k.saturating_add(dk).min(nz-1),
        k.saturating_sub(dk).max(0),
    );
    let f1 = gamma[ip][jp][kp][a][b];
    let f2 = gamma[im][jm][km][a][b];
    (f1 - f2) / (2.0 * h)
}

/// Helper: compute spatial derivative of Christoffel symbol at (i,j,k) with respect to dir (0=x,1=y,2=z)
fn spatial_deriv_christoffel(
    christoffel: &[[[f64; 3]; 3]; 3],
    i: usize, j: usize, k: usize,
    dir: usize, a: usize, b: usize, c: usize,
    dx: f64, dy: f64, dz: f64,
    nx: usize, ny: usize, nz: usize,
) -> f64 {
    // For now, use the same point (no grid for Christoffel, so just return 0)
    // In a full implementation, Christoffel symbols would be stored on the grid and this would use finite differences
    0.0
}
/// Utility: finite difference for first spatial derivative (central difference, periodic BC)
fn fd1(field: &Vec<f64>, n: usize, stride: usize, i: usize, dx: f64) -> f64 {
    // stride: 1 for x, nx for y, nx*ny for z
    let ip = if i + stride < n { i + stride } else { i + stride - n };
    let im = if i >= stride { i - stride } else { n + i - stride };
    (field[ip] - field[im]) / (2.0 * dx)
}

/// Utility: finite difference for second spatial derivative (central difference, periodic BC)
fn fd2(field: &Vec<f64>, n: usize, stride: usize, i: usize, dx: f64) -> f64 {
    let ip = if i + stride < n { i + stride } else { i + stride - n };
    let im = if i >= stride { i - stride } else { n + i - stride };
    (field[ip] - 2.0 * field[i] + field[im]) / (dx * dx)
}
// numerical relativity infrastructure
//
// this module provides the data structures and stubs for evolving the spacetime metric numerically on a grid.
// it is the foundation for a full numerical relativity solver (e.g., BSSN or ADM formalism).

use crate::relativity::FourVector;

/// 3+1D grid for metric field and auxiliary variables
#[derive(Clone)]
pub struct Grid3D {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
}

impl Grid3D {
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> Self {
        Self { nx, ny, nz, dx, dy, dz }
    }
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nx * (j + self.ny * k)
    }
}

/// symmetric 4x4 metric tensor at each grid point
type MetricTensor = [[f64; 4]; 4];


/// 3-metric (spatial part of metric tensor)
pub type Metric3 = [[f64; 3]; 3];
/// extrinsic curvature tensor (symmetric 3x3)
pub type ExtrinsicCurvature = [[f64; 3]; 3];


/// Main struct for numerical relativity gravity (metric evolution)
// Derive Clone for NumericalRelativityGravity so it can be cloned for model comparison
#[derive(Clone)]
pub struct NumericalRelativityGravity {
    pub grid: Grid3D,
    pub metric: Vec<Metric3>, // 3-metric at each grid point
    pub extrinsic_curvature: Vec<ExtrinsicCurvature>,
    pub lapse: Vec<f64>, // lapse function at each grid point
    pub shift: Vec<[f64; 3]>, // shift vector at each grid point
    pub shift_aux: Vec<[f64; 3]>, // B^i auxiliary field for Gamma-driver
    // --- BSSN variables ---
    pub conformal_metric: Vec<Metric3>, // \bar{\gamma}_{ij}
    pub conformal_factor: Vec<f64>,     // \phi (logarithmic, so \gamma_{ij} = e^{4\phi} \bar{\gamma}_{ij})
    pub trace_free_extrinsic: Vec<ExtrinsicCurvature>, // \bar{A}_{ij}
    pub trace_K: Vec<f64>, // K (trace of extrinsic curvature)
    pub conformal_connection: Vec<[f64; 3]>, // \bar{\Gamma}^i
}

impl NumericalRelativityGravity {
    pub fn new(grid: Grid3D) -> Self {
        let n = grid.nx * grid.ny * grid.nz;
        // initialize to flat space (Euclidean 3-metric, zero extrinsic curvature, lapse=1, shift=0)
        let flat3 = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let zero_k = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]];
        Self {
            grid,
            metric: vec![flat3; n],
            extrinsic_curvature: vec![zero_k; n],
            lapse: vec![1.0; n],
            shift: vec![[0.0, 0.0, 0.0]; n],
            shift_aux: vec![[0.0, 0.0, 0.0]; n],
            // --- BSSN variables ---
            conformal_metric: vec![flat3; n],
            conformal_factor: vec![0.0; n],
            trace_free_extrinsic: vec![zero_k; n],
            trace_K: vec![0.0; n],
            conformal_connection: vec![[0.0; 3]; n],
        }
    }

    /// Set initial data for a grid point
    pub fn set_initial_data(&mut self, i: usize, j: usize, k: usize, metric: Metric3, k_tensor: ExtrinsicCurvature, lapse: f64, shift: [f64; 3]) {
        let idx = self.grid.index(i, j, k);
        self.metric[idx] = metric;
        self.extrinsic_curvature[idx] = k_tensor;
        self.lapse[idx] = lapse;
        self.shift[idx] = shift;
    }

    /// Get the 3-metric at a given spatial index
    pub fn metric_at(&self, i: usize, j: usize, k: usize) -> &Metric3 {
        &self.metric[self.grid.index(i, j, k)]
    }

    /// Get the extrinsic curvature at a given spatial index
    pub fn extrinsic_curvature_at(&self, i: usize, j: usize, k: usize) -> &ExtrinsicCurvature {
        &self.extrinsic_curvature[self.grid.index(i, j, k)]
    }

    /// Get the lapse at a given spatial index
    pub fn lapse_at(&self, i: usize, j: usize, k: usize) -> f64 {
        self.lapse[self.grid.index(i, j, k)]
    }

    /// Get the shift at a given spatial index
    pub fn shift_at(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        self.shift[self.grid.index(i, j, k)]
    }

    /// Evolve the metric, extrinsic curvature, lapse, and shift by one time step (ADM update)
    pub fn evolve_step(&mut self, dt: f64) {
        let eta = 1.0; // Damping parameter for shift evolution (can be made configurable)
        let gamma_driver_eta = 2.0; // Gamma-driver parameter
        let lapse_damping = 0.2; // Constraint damping for lapse
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let nz = self.grid.nz;
        let n = nx * ny * nz;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        let dz = self.grid.dz;

        // Evolve the lapse using Bona-Masso slicing (1+log slicing) with constraint damping
        // ∂_t α = -2 α K + ζ Δα
        let mut new_lapse = self.lapse.clone();
        // Prepare 3D array for Laplacian
        let mut lapse_3d = vec![vec![vec![0.0; nz]; ny]; nx];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.grid.index(i, j, k);
                    lapse_3d[i][j][k] = self.lapse[idx];
                }
            }
        }
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.grid.index(i, j, k);
                    let K = self.extrinsic_curvature[idx][0][0] + self.extrinsic_curvature[idx][1][1] + self.extrinsic_curvature[idx][2][2];
                    // Laplacian Δα ≈ (α_{i+1,j,k} + α_{i-1,j,k} - 2α_{i,j,k})/dx^2 + ...
                    let lap = {
                        let ip = if i+1 < nx { i+1 } else { 0 };
                        let im = if i > 0 { i-1 } else { nx-1 };
                        let jp = if j+1 < ny { j+1 } else { 0 };
                        let jm = if j > 0 { j-1 } else { ny-1 };
                        let kp = if k+1 < nz { k+1 } else { 0 };
                        let km = if k > 0 { k-1 } else { nz-1 };
                        (lapse_3d[ip][j][k] + lapse_3d[im][j][k] - 2.0*lapse_3d[i][j][k])/(dx*dx)
                        + (lapse_3d[i][jp][k] + lapse_3d[i][jm][k] - 2.0*lapse_3d[i][j][k])/(dy*dy)
                        + (lapse_3d[i][j][kp] + lapse_3d[i][j][km] - 2.0*lapse_3d[i][j][k])/(dz*dz)
                    };
                    new_lapse[idx] += (-2.0 * self.lapse[idx] * K + lapse_damping * lap) * dt;
                }
            }
        }
        self.lapse = new_lapse;

        // --- Gamma-driver shift evolution ---
        // Compute conformal connection functions (approximate, using metric derivatives)
        // For now, use finite differences on the metric diagonal as a proxy for Γ^i
        let mut gamma_conn = vec![[0.0; 3]; n];
        for k in 1..(nz-1) {
            for j in 1..(ny-1) {
                for i in 1..(nx-1) {
                    let idx = self.grid.index(i, j, k);
                    // Γ^i ≈ sum_j ∂_j γ^{ij} (very rough, for demonstration)
                    for d in 0..3 {
                        let mut sum = 0.0;
                        // Central difference for each direction
                        if d == 0 {
                            let ip = self.grid.index(i+1, j, k);
                            let im = self.grid.index(i-1, j, k);
                            sum += (self.metric[ip][d][d] - self.metric[im][d][d]) / (2.0 * dx);
                        } else if d == 1 {
                            let jp = self.grid.index(i, j+1, k);
                            let jm = self.grid.index(i, j-1, k);
                            sum += (self.metric[jp][d][d] - self.metric[jm][d][d]) / (2.0 * dy);
                        } else if d == 2 {
                            let kp = self.grid.index(i, j, k+1);
                            let km = self.grid.index(i, j, k-1);
                            sum += (self.metric[kp][d][d] - self.metric[km][d][d]) / (2.0 * dz);
                        }
                        gamma_conn[idx][d] = sum;
                    }
                }
            }
        }
        // Evolve B^i and shift β^i: ∂_t β^i = (3/4) B^i, ∂_t B^i = ∂_t Γ^i - η B^i
        let mut new_shift = self.shift.clone();
        let mut new_shift_aux = self.shift_aux.clone();
        for idx in 0..n {
            for d in 0..3 {
                // For demonstration, use gamma_conn as ∂_t Γ^i (should be time derivative, but use spatial for now)
                let dgamma = gamma_conn[idx][d];
                new_shift_aux[idx][d] += (dgamma - gamma_driver_eta * self.shift_aux[idx][d]) * dt;
                new_shift[idx][d] += 0.75 * new_shift_aux[idx][d] * dt;
            }
        }
        self.shift = new_shift;
        self.shift_aux = new_shift_aux;

        // Gather 3-metric, extrinsic curvature, and lapse into 4D arrays for evolution
        let mut gamma = vec![vec![vec![[[0.0; 3]; 3]; nz]; ny]; nx];
        let mut K = vec![vec![vec![[[0.0; 3]; 3]; nz]; ny]; nx];
        let mut lapse = vec![vec![vec![0.0; nz]; ny]; nx];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.grid.index(i, j, k);
                    gamma[i][j][k] = self.metric[idx];
                    K[i][j][k] = self.extrinsic_curvature[idx];
                    lapse[i][j][k] = self.lapse[idx];
                }
            }
        }

        // Evolve metric and extrinsic curvature using ADM update
        evolve_adm_euler(&mut gamma, &mut K, &lapse, dx, dy, dz, dt);

        // Apply periodic boundary conditions (can be made runtime-configurable)
        apply_periodic_boundary_conditions(&mut gamma, &mut K);

        // Write updated values back to flat storage
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.grid.index(i, j, k);
                    self.metric[idx] = gamma[i][j][k];
                    self.extrinsic_curvature[idx] = K[i][j][k];
                }
            }
        }

        // --- Update BSSN variables from ADM fields ---
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.grid.index(i, j, k);
                    let g = &self.metric[idx];
                    // Compute conformal factor phi = (1/12) * ln(det(gamma_ij))
                    let det = g[0][0]*g[1][1]*g[2][2] + 2.0*g[0][1]*g[0][2]*g[1][2]
                        - g[0][0]*g[1][2]*g[1][2] - g[1][1]*g[0][2]*g[0][2] - g[2][2]*g[0][1]*g[0][1];
                    let phi = (det.abs().ln() / 12.0).max(-20.0).min(20.0); // Clamp for safety
                    self.conformal_factor[idx] = phi;
                    // Compute conformal metric \bar{\gamma}_{ij} = e^{-4phi} gamma_{ij}
                    let exp_m4phi = (-4.0*phi).exp();
                    let mut conf_g = [[0.0; 3]; 3];
                    for a in 0..3 {
                        for b in 0..3 {
                            conf_g[a][b] = g[a][b] * exp_m4phi;
                        }
                    }
                    self.conformal_metric[idx] = conf_g;
                    // Compute trace K
                    let Kij = &self.extrinsic_curvature[idx];
                    let K_trace = Kij[0][0] + Kij[1][1] + Kij[2][2];
                    self.trace_K[idx] = K_trace;
                    // Compute trace-free extrinsic curvature \bar{A}_{ij} = e^{-4phi}(K_{ij} - 1/3 gamma_{ij} K)
                    let mut Abar = [[0.0; 3]; 3];
                    for a in 0..3 {
                        for b in 0..3 {
                            Abar[a][b] = exp_m4phi * (Kij[a][b] - (1.0/3.0)*g[a][b]*K_trace);
                        }
                    }
                    self.trace_free_extrinsic[idx] = Abar;
                    // Compute conformal connection functions \bar{\Gamma}^i = -\partial_j \bar{\gamma}^{ij}
                    // For now, use central differences on the diagonal (approximate)
                    let mut inv_conf_g = invert_3x3(&conf_g).unwrap_or([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]);
                    let mut gamma_bar = [0.0; 3];
                    for d in 0..3 {
                        let mut sum = 0.0;
                        if d == 0 && i > 0 && i+1 < nx {
                            let idxp = self.grid.index(i+1, j, k);
                            let idxm = self.grid.index(i-1, j, k);
                            sum += (inv_conf_g[d][d] - invert_3x3(&self.conformal_metric[idxm]).unwrap_or(inv_conf_g)[d][d]) / (2.0*dx);
                        } else if d == 1 && j > 0 && j+1 < ny {
                            let idxp = self.grid.index(i, j+1, k);
                            let idxm = self.grid.index(i, j-1, k);
                            sum += (inv_conf_g[d][d] - invert_3x3(&self.conformal_metric[idxm]).unwrap_or(inv_conf_g)[d][d]) / (2.0*dy);
                        } else if d == 2 && k > 0 && k+1 < nz {
                            let idxp = self.grid.index(i, j, k+1);
                            let idxm = self.grid.index(i, j, k-1);
                            sum += (inv_conf_g[d][d] - invert_3x3(&self.conformal_metric[idxm]).unwrap_or(inv_conf_g)[d][d]) / (2.0*dz);
                        }
                        gamma_bar[d] = -sum;
                    }
                    self.conformal_connection[idx] = gamma_bar;
                }
            }
        }

        // --- Constraint monitoring diagnostics ---
        let h_constraint = compute_hamiltonian_constraint(&gamma, &K, dx, dy, dz);
        let mut h_max: f64 = 0.0;
        let mut h_sum: f64 = 0.0;
        let mut h_count: f64 = 0.0;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let val = h_constraint[i][j][k].abs();
                    h_max = h_max.max(val);
                    h_sum += val;
                    h_count += 1.0;
                }
            }
        }
        let h_avg = h_sum / h_count.max(1.0);
        println!("[ADM evolve] Hamiltonian constraint: max={:.3e}, avg={:.3e}", h_max, h_avg);
        // TODO: Add momentum constraint diagnostics
        // NOTE: This is a minimal ADM update. Real evolution requires constraint handling and boundary conditions.
    }
}
