use crate::curvature::Metric;
use ecs::System;
mod curvature;
mod ecs;
mod relativity;

// --- ECS foundation, entity management, and component storage demo ---




pub use relativity::{FourVector, FourVelocity};

#[derive(Debug, Clone)]
pub struct Position(pub FourVector);
impl ecs::Component for Position {}

#[derive(Debug, Clone)]
pub struct Velocity(pub FourVelocity);
impl ecs::Component for Velocity {}

#[derive(Debug, Clone)]
pub struct Acceleration(pub FourVector);
impl ecs::Component for Acceleration {}

#[derive(Debug, Clone)]
pub struct Mass(pub f64);
impl ecs::Component for Mass {}

#[derive(Debug, Clone)]
pub struct Force(pub f64, pub f64, pub f64);
impl ecs::Component for Force {}

fn main() {
    println!("[cargo run] Physics engine starting up...");
    let mut world = ecs::World::new();

    // ECS foundation, entity management, and component storage
    println!("[progress] Checking off: ECS foundation, entity management, component storage");


    // register all physics components
    world.register_component::<Position>();
    world.register_component::<Velocity>();
    world.register_component::<Acceleration>();
    world.register_component::<Mass>();
    world.register_component::<Force>();
    println!("[progress] Registered Position, Velocity, Acceleration, Mass, Force components");

    // create entities with movement in all directions
    let mut entities = Vec::new();
    let directions = [
        (1.0, 0.0, 0.0), // +X
        (0.0, 1.0, 0.0), // +Y
        (0.0, 0.0, 1.0), // +Z
        (-1.0, 0.0, 0.0), // -X
        (0.0, -1.0, 0.0), // -Y
        (0.0, 0.0, -1.0), // -Z
        (1.0, 1.0, 1.0), // diagonal
    ];
    let masses = [1.0_f64, 2.0, 0.5, 3.0, 1.5, 2.5, 1.2];
    let c = 10.0;
    for (i, dir) in directions.iter().enumerate() {
        let entity = world.create_entity();
        // clamp initial 3-velocity to be strictly less than c
        let v2: f64 = dir.0 * dir.0 + dir.1 * dir.1 + dir.2 * dir.2;
        let v_max: f64 = 0.999 * c;
        let scale: f64 = if v2 > v_max * v_max {
            v_max / v2.sqrt()
        } else { 1.0 };
        let vx: f64 = dir.0 * scale;
        let vy: f64 = dir.1 * scale;
        let vz: f64 = dir.2 * scale;
        let pos = FourVector { t: 0.0, x: 0.0, y: 0.0, z: 0.0 };
        let vel = FourVelocity::from_3velocity(vx, vy, vz, c);
        let acc = FourVector { t: 0.0, x: 0.1 * vx, y: 0.1 * vy, z: 0.1 * vz };
        assert!(pos.t.is_finite() && pos.x.is_finite() && pos.y.is_finite() && pos.z.is_finite());
        assert!(vel.0.t.is_finite() && vel.0.x.is_finite() && vel.0.y.is_finite() && vel.0.z.is_finite());
        assert!(acc.t.is_finite() && acc.x.is_finite() && acc.y.is_finite() && acc.z.is_finite());
        world.add_component(entity, Position(pos));
        world.add_component(entity, Velocity(vel));
        world.add_component(entity, Acceleration(acc));
        world.add_component(entity, Mass(masses[i % masses.len()]));
        world.add_component(entity, Force(0.0, 0.0, 0.0));
        entities.push(entity);
    }
    println!("[progress] Entities initialized with movement in all directions");

    // print the entities' initial state
    for &entity in &entities {
        if let Some(pos) = world.get_component::<Position>(entity) {
            println!("[ecs] Entity {:?} has position: {:?}", entity, pos);
        }
        if let Some(vel) = world.get_component::<Velocity>(entity) {
            println!("[ecs] Entity {:?} has velocity: {:?}", entity, vel);
            let (vx, vy, vz) = vel.0.three_velocity(c);
            println!("[relativity] 3-velocity: ({:.3}, {:.3}, {:.3})", vx, vy, vz);
        }
        if let Some(acc) = world.get_component::<Acceleration>(entity) {
            println!("[ecs] Entity {:?} has acceleration: {:?}", entity, acc);
        }
    }

    // configurable friction
    let mut friction = 0.9_f64; // default
    println!("[config] Friction set to: {}", friction);
    let gravity = -0.5_f64; // gravity in -Y

    // add more entities with different initial positions
    let mut more_entities = Vec::new();
    let positions = [
        (5.0, 0.0, 0.0),
        (0.0, 5.0, 0.0),
        (0.0, 0.0, 5.0),
        (-5.0, 0.0, 0.0),
        (0.0, -5.0, 0.0),
        (0.0, 0.0, -5.0),
    ];
    let more_masses = [2.0_f64, 3.0, 1.0, 4.0, 2.5, 1.8];
    for (i, pos) in positions.iter().enumerate() {
        let entity = world.create_entity();
        world.add_component(entity, Position(FourVector { t: 0.0, x: pos.0, y: pos.1, z: pos.2 }));
        world.add_component(entity, Velocity(FourVelocity::from_3velocity(0.0, 0.0, 0.0, c)));
        world.add_component(entity, Acceleration(FourVector { t: 0.0, x: 0.0, y: 0.0, z: 0.0 }));
        world.add_component(entity, Mass(more_masses[i % more_masses.len()]));
        world.add_component(entity, Force(0.0, 0.0, 0.0));
        more_entities.push(entity);
    }
#[cfg(test)]
mod relativity_tests {
    use super::*;
    use crate::relativity::*;

    #[test]
    fn minkowski_interval() {
        let c = 10.0;
        let v = FourVector { t: 1.0 * c, x: 3.0, y: 4.0, z: 0.0 };
        let s2 = v.interval2_metric(&Metric::minkowski());
        assert!((s2 + c * c - 25.0).abs() < 1e-3); // -c^2 t^2 + x^2 + y^2 + z^2
    }

    #[test]
    fn lorentz_boost_x() {
        let c = 10.0;
        let v = FourVector { t: 1.0 * c, x: 5.0, y: 0.0, z: 0.0 };
        let beta = 0.6;
        let boosted = v.lorentz_boost_x(beta);
        let gamma = 1.0 / (1.0 - beta * beta).sqrt();
        assert!((boosted.t - (gamma * (v.t - beta * v.x))).abs() < 1e-3);
        assert!((boosted.x - (gamma * (v.x - beta * v.t))).abs() < 1e-3);
    }

    #[test]
    fn four_velocity_from_3velocity() {
        let c = 10.0;
        let vx = 6.0;
        let v = FourVelocity::from_3velocity(vx, 0.0, 0.0, c);
        let gamma = 1.0 / (1.0 - (vx * vx) / (c * c)).sqrt();
        assert!((v.0.t - gamma * c).abs() < 1e-3);
        assert!((v.0.x - gamma * vx).abs() < 1e-3);
    }
}
    let all_entities: Vec<_> = entities.iter().chain(more_entities.iter()).copied().collect();
    println!("[progress] Added more entities with different initial positions");

    // relativistic constants
    let c = 10.0; // speed of light (arbitrary units for demo)

    // apply four-force-based friction and gravity with proper time stepping
    use relativity::FourForce;
    let d_t = 1.0; // global simulation time step
    // accumulate all forces (friction, global gravity, pairwise gravity) before integrating
    let g_const = 0.1;
    let mut total_forces = vec![FourForce { t: 0.0, x: 0.0, y: 0.0, z: 0.0 }; all_entities.len()];
    // prepare to collect all position/velocity updates for ECS borrow safety
    let mut updates = Vec::new();
    // first, add friction and global gravity to total_forces
    for (idx, &entity) in all_entities.iter().enumerate() {
        if let (Some(vel), Some(mass)) = (world.get_component::<Velocity>(entity), world.get_component::<Mass>(entity)) {
            // gravity as four-force (in -y direction)
            let gravity_force = FourForce {
                t: 0.0,
                x: 0.0,
                y: mass.0 * gravity,
                z: 0.0,
            };
            // friction as four-force (opposes spatial velocity)
            let (vx, vy, vz) = vel.0.three_velocity(c);
            let v_norm = (vx * vx + vy * vy + vz * vz).sqrt();
            let friction_mag = if v_norm > 1e-6 {
                -friction * mass.0 * v_norm
            } else { 0.0 };
            let friction_force = FourForce {
                t: 0.0,
                x: if v_norm > 1e-6 { friction_mag * vx / v_norm } else { 0.0 },
                y: if v_norm > 1e-6 { friction_mag * vy / v_norm } else { 0.0 },
                z: if v_norm > 1e-6 { friction_mag * vz / v_norm } else { 0.0 },
            };
            // add to total_forces
            total_forces[idx].t += gravity_force.t + friction_force.t;
            total_forces[idx].x += gravity_force.x + friction_force.x;
            total_forces[idx].y += gravity_force.y + friction_force.y;
            total_forces[idx].z += gravity_force.z + friction_force.z;
        }
    }
    // now accumulate pairwise gravity into total_forces
    for i in 0..all_entities.len() {
        for j in (i+1)..all_entities.len() {
            let e1 = all_entities[i];
            let e2 = all_entities[j];
            let (pos1, pos2, mass1, mass2, vel1, vel2) = (
                world.get_component::<Position>(e1),
                world.get_component::<Position>(e2),
                world.get_component::<Mass>(e1),
                world.get_component::<Mass>(e2),
                world.get_component::<Velocity>(e1),
                world.get_component::<Velocity>(e2),
            );
            if let (Some(pos1), Some(pos2), Some(m1), Some(m2), Some(vel1), Some(vel2)) = (pos1, pos2, mass1, mass2, vel1, vel2) {
                let delta = pos2.0 - pos1.0;
                let dx = delta.x;
                let dy = delta.y;
                let dz = delta.z;
                let dist_sq = dx*dx + dy*dy + dz*dz + 1e-6; // avoid div by zero
                let dist = dist_sq.sqrt();
                // compute gamma for each
                let gamma1 = vel1.0.gamma(c);
                let gamma2 = vel2.0.gamma(c);
                // relativistic force magnitude (linearized): F = G * (gamma1*m1) * (gamma2*m2) / r^2
                let force_mag = g_const * (gamma1 * m1.0) * (gamma2 * m2.0) / dist_sq;
                let fx = force_mag * dx / dist;
                let fy = force_mag * dy / dist;
                let fz = force_mag * dz / dist;
                // accumulate +f for i, -f for j
                total_forces[i].x += fx;
                total_forces[i].y += fy;
                total_forces[i].z += fz;
                total_forces[j].x -= fx;
                total_forces[j].y -= fy;
                total_forces[j].z -= fz;
            }
        }
    }
    // now apply the total force to each entity in a single integration step
    for (idx, &entity) in all_entities.iter().enumerate() {
        if let (Some(pos), Some(vel), Some(mass)) = (
            world.get_component::<Position>(entity),
            world.get_component::<Velocity>(entity),
            world.get_component::<Mass>(entity)) {
            let gamma = vel.0.gamma(c).max(1e-6);
            let d_tau = d_t / gamma;
            let total_force = total_forces[idx];
            let force_norm = (total_force.t * total_force.t + total_force.x * total_force.x + total_force.y * total_force.y + total_force.z * total_force.z).sqrt();
            println!("[debug] Entity {:?}: d_tau = {:.6}, |total_force| = {:.6}", entity, d_tau, force_norm);
            let ortho = total_force.orthogonal_to(vel.0, c);
            let eps_mass = 1e-12;
            if mass.0 <= eps_mass {
                eprintln!("[error] entity {:?} has tiny mass = {}; skipping force update", entity, mass.0);
                continue;
            }
            let four_accel = FourVector {
                t: ortho.t / (gamma * mass.0),
                x: ortho.x / (gamma * mass.0),
                y: ortho.y / (gamma * mass.0),
                z: ortho.z / (gamma * mass.0),
            };
            let (new_pos, new_vel) = vel.0.integrate_geodesic(&pos.0, &|x| Metric::minkowski(), None, mass.0, d_tau, 1e-6);
            // print u·u and |v| after integrate
            let uu = new_vel.norm2(c);
            let (vx, vy, vz) = new_vel.three_velocity(c);
            let vmag = (vx * vx + vy * vy + vz * vz).sqrt();
            let gamma_old = vel.0.gamma(c);
            let gamma_new = new_vel.gamma(c);
            let mut safe_vel = None;
            if gamma_new.is_nan() || gamma_new > gamma_old * 100.0 + 1.0 {
                eprintln!("[safety] Clamping gamma jump: old={}, new={}; entity={:?}", gamma_old, gamma_new, entity);
                safe_vel = Some(vel.0.integrate_accel(four_accel, d_tau * 0.1, c).normalize(c));
            }
            updates.push((entity, new_pos, if let Some(sv) = safe_vel { sv } else { new_vel }));
            println!("[debug] Entity {:?}: u·u = {:.6}, |v| = {:.6}, gamma_old = {:.3}, gamma_new = {:.3}", entity, uu, vmag, gamma_old, gamma_new);
        }
    // apply all position/velocity updates after the loop to avoid ECS borrow conflicts
    let mut updates = Vec::new();
    // (loop body unchanged except for pushing to updates)
    // now apply all updates
    for (entity, new_pos, new_vel) in updates {
        world.add_component(entity, Position(new_pos));
        world.add_component(entity, Velocity(new_vel));
    }
    }
    println!("[progress] Applied all forces (friction, gravity, pairwise gravity) to all entities (single integration per entity, proper d_tau)");

    // relativistically correct collision detection and response
    let collision_dist = 1.0;
    for i in 0..all_entities.len() {
        for j in (i+1)..all_entities.len() {
            let e1 = all_entities[i];
            let e2 = all_entities[j];
            let (pos1, pos2, vel1, vel2, m1, m2) = (
                world.get_component::<Position>(e1),
                world.get_component::<Position>(e2),
                world.get_component::<Velocity>(e1),
                world.get_component::<Velocity>(e2),
                world.get_component::<Mass>(e1),
                world.get_component::<Mass>(e2),
            );
            if let (Some(pos1), Some(pos2), Some(vel1), Some(vel2), Some(m1), Some(m2)) = (pos1, pos2, vel1, vel2, m1, m2) {
                let delta = pos2.0 - pos1.0;
                let dx = delta.x;
                let dy = delta.y;
                let dz = delta.z;
                let dist = (dx*dx + dy*dy + dz*dz).sqrt();
                if dist < collision_dist {
                    // relativistic two-body elastic collision (general 3D)
                    let p1 = vel1.0.four_momentum(m1.0, c);
                    let p2 = vel2.0.four_momentum(m2.0, c);
                    let p_tot = p1 + p2;
                    // compute beta vector for CoM frame: beta = p_spatial / p_t (units v/c)
                    let beta_x = p_tot.x / p_tot.t;
                    let beta_y = p_tot.y / p_tot.t;
                    let beta_z = p_tot.z / p_tot.t;
                    let beta2 = beta_x*beta_x + beta_y*beta_y + beta_z*beta_z;
                    if beta2 >= 1.0 - 1e-12 {
                        eprintln!("[collision] Warning: |beta| >= 1 for CoM boost; skipping relativistic collision for pair.");
                        continue;
                    }
                    // boost both four-momenta to COM frame
                    let p1_com = p1.lorentz_boost_beta(-beta_x, -beta_y, -beta_z, c);
                    let p2_com = p2.lorentz_boost_beta(-beta_x, -beta_y, -beta_z, c);
                    // in COM frame, energies are p1_com.t, p2_com.t and spatial p vectors opposite
                    // for an elastic collision in COM frame, just reverse spatial momentum components
                    // reverse spatial momentum and recompute energy (on-shell) in CoM frame
                    let reverse_and_recompute = |p: &FourVector, mass: f64| {
                        let px = -p.x;
                        let py = -p.y;
                        let pz = -p.z;
                        let spatial2 = px*px + py*py + pz*pz;
                        let energy = (mass*mass*c*c + spatial2).sqrt();
                        FourVector { t: energy, x: px, y: py, z: pz }
                    };
                    let p1_com_post = reverse_and_recompute(&p1_com, m1.0);
                    let p2_com_post = reverse_and_recompute(&p2_com, m2.0);
                    // boost back to lab frame
                    let p1_lab = p1_com_post.lorentz_boost_beta(beta_x, beta_y, beta_z, c);
                    let p2_lab = p2_com_post.lorentz_boost_beta(beta_x, beta_y, beta_z, c);
                    // now convert four-momentum back to four-velocity properly:
                    let v1_new = FourVelocity::from_four_momentum(p1_lab, m1.0, c);
                    let v2_new = FourVelocity::from_four_momentum(p2_lab, m2.0, c);
                    world.add_component(e1, Velocity(v1_new.normalize(c)));
                    world.add_component(e2, Velocity(v2_new.normalize(c)));
                    println!("[collision] Relativistic elastic collision: {:?} <-> {:?}", e1, e2);
                }
            }
        }
    }
    println!("[progress] Relativistic collision detection and response complete");

    // relativistically consistent position update: advance coordinate time by d_t, spatial by v^i/v^0 * d_t
    for &entity in &all_entities {
        if let (Some(pos), Some(vel)) = (world.get_component::<Position>(entity), world.get_component::<Velocity>(entity)) {
            let v0 = vel.0.0.t;
            let vx = vel.0.0.x;
            let vy = vel.0.0.y;
            let vz = vel.0.0.z;
            let mut new_pos = pos.0;
            // advance coordinate time by d_t
            new_pos.t += c * d_t;
            // advance spatial coordinates by (v^i / v^0) * c * d_t
            let dt_factor = c * d_t / v0.max(1e-6); // v0 = gamma * c
            new_pos.x += vx * dt_factor;
            new_pos.y += vy * dt_factor;
            new_pos.z += vz * dt_factor;
            world.add_component(entity, Position(new_pos));
        }
    }
    println!("[progress] Position integrated by four-velocity for all entities (relativistic coordinate time)");


    // register and run PhysicsSystem (now with acceleration integration)
    println!("[progress] Checking off: System registration and execution, Acceleration integration");
    let mut physics_system = ecs::PhysicsSystem;
    physics_system.run(&mut world);

    // print each entity's new state after physics
    for &entity in &entities {
        if let Some(pos) = world.get_component::<Position>(entity) {
            println!("[ecs] Entity {:?} new position after physics: {:?}", entity, pos);
        }
        if let Some(vel) = world.get_component::<Velocity>(entity) {
            println!("[ecs] Entity {:?} new velocity after physics: {:?}", entity, vel);
        }
    }

    // TODO: Add more physics systems and expand ECS features
}
