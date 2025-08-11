/// example system: updates positions by velocity
pub struct PhysicsSystem;

impl System for PhysicsSystem {
    fn run(&mut self, world: &mut World) {
        println!("[system] Running PhysicsSystem (integrates Acceleration -> Velocity -> Position)");
        let type_id_pos = std::any::TypeId::of::<crate::Position>();
        let type_id_vel = std::any::TypeId::of::<crate::Velocity>();
        let type_id_acc = std::any::TypeId::of::<crate::Acceleration>();

        // 1. integrate acceleration into velocity
        // 1. integrate four-acceleration into four-velocity (proper time step = 1.0 for demo)
        let acc_vec: Vec<(Entity, crate::Acceleration)> =
            if let Some(acc_map) = world.components.get(&type_id_acc)
                .and_then(|c| c.downcast_ref::<std::collections::HashMap<Entity, crate::Acceleration>>()) {
                acc_map.iter().map(|(&e, a)| (e, a.clone())).collect()
            } else { vec![] };
        if let Some(vel_map) = world.components.get_mut(&type_id_vel)
            .and_then(|c| c.downcast_mut::<std::collections::HashMap<Entity, crate::Velocity>>()) {
            let c = 10.0; // use same c as main.rs for demo
            let d_t = 1.0; // use same d_t as main.rs for demo
            for (entity, acc) in acc_vec {
                if let Some(vel) = vel_map.get_mut(&entity) {
                    let gamma = vel.0.gamma(c).max(1e-6);
                    let d_tau = d_t / gamma;
                    // when integrating four-acceleration, use geodesic integrator if in curved spacetime:
                    // let (new_pos, new_vel) = vel.0.integrate_geodesic(&pos.0, &|x| Metric::minkowski(), None, mass.0, d_tau, 1e-6);
                    vel.0 = vel.0.integrate_accel(acc.0, d_tau, c).normalize(c);
                    println!("[physics] Entity {:?} velocity updated by four-acceleration: {:?}", entity, vel);
                }
            }
        }
        println!("[progress] Four-acceleration integrated into four-velocity");

        // 2. integrate four-velocity into position (proper time step = 1.0 for demo)
        let vel_vec: Vec<(Entity, crate::Velocity)> =
            if let Some(vel_map) = world.components.get(&type_id_vel)
                .and_then(|c| c.downcast_ref::<std::collections::HashMap<Entity, crate::Velocity>>()) {
                vel_map.iter().map(|(&e, v)| (e, v.clone())).collect()
            } else { vec![] };
        if let Some(pos_map) = world.components.get_mut(&type_id_pos)
            .and_then(|c| c.downcast_mut::<std::collections::HashMap<Entity, crate::Position>>()) {
            for (entity, vel) in vel_vec {
                if let Some(pos) = pos_map.get_mut(&entity) {
                    // for all four-vector math, use metric where appropriate.
                    pos.0.t += vel.0.0.t * 1.0;
                    pos.0.x += vel.0.0.x * 1.0;
                    pos.0.y += vel.0.0.y * 1.0;
                    pos.0.z += vel.0.0.z * 1.0;
                    println!("[physics] Entity {:?} moved to position: {:?}", entity, pos);
                }
            }
        }
        println!("[progress] Four-velocity integrated into position");
    }
}
/// ECS core traits and types for the engine

use std::collections::{HashMap, HashSet};
use std::any::{Any, TypeId};

/// unique identifier for an entity.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Entity(pub u32);

/// trait for all components.
pub trait Component: Send + Sync {}

/// trait for all systems.
pub trait System {
    fn run(&mut self, world: &mut World);
}

/// the world struct holds all entities and components.
pub struct World {
    next_entity: u32,
    alive_entities: HashSet<Entity>,
    components: HashMap<TypeId, Box<dyn Any + Send + Sync>>,
}

impl World {
    /// create a new, empty world
    pub fn new() -> Self {
        println!("[ecs] World created (entity management, component storage initialized)");
        Self {
            next_entity: 0,
            alive_entities: HashSet::new(),
            components: HashMap::new(),
        }
    }

    /// create a new entity and return its id
    pub fn create_entity(&mut self) -> Entity {
        let id = self.next_entity;
        self.next_entity += 1;
        let entity = Entity(id);
        self.alive_entities.insert(entity);
        println!("[ecs] Entity {:?} created", entity);
        entity
    }

    /// delete an entity and remove its components
    pub fn delete_entity(&mut self, entity: Entity) {
        self.alive_entities.remove(&entity);
        // TODO: Remove components for this entity
        println!("[ecs] Entity {:?} deleted", entity);
    }

    /// register a component type
    pub fn register_component<T: Component + 'static>(&mut self) {
        let type_id = TypeId::of::<T>();
        if !self.components.contains_key(&type_id) {
            self.components.insert(type_id, Box::new(HashMap::<Entity, T>::new()));
            println!("[ecs] Registered component type: {}", std::any::type_name::<T>());
        }
    }

    /// add a component to an entity
    pub fn add_component<T: Component + 'static>(&mut self, entity: Entity, component: T) {
        let type_id = TypeId::of::<T>();
        if let Some(storage) = self.components.get_mut(&type_id) {
            let storage = storage.downcast_mut::<HashMap<Entity, T>>().unwrap();
            let existed = storage.contains_key(&entity);
            storage.insert(entity, component);
            if existed {
                println!("[ecs] Replaced component {} for entity {:?}", std::any::type_name::<T>(), entity);
            } else {
                println!("[ecs] Added component {} to entity {:?}", std::any::type_name::<T>(), entity);
            }
        } else {
            panic!("Component type not registered");
        }
    }

    /// get a component for an entity
    pub fn get_component<T: Component + 'static>(&self, entity: Entity) -> Option<&T> {
        let type_id = TypeId::of::<T>();
        self.components.get(&type_id)
            .and_then(|storage| storage.downcast_ref::<HashMap<Entity, T>>())
            .and_then(|map| map.get(&entity))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entity_id_unique() {
        let e1 = Entity(1);
        let e2 = Entity(2);
        assert_ne!(e1, e2);
    }

    #[test]
    fn world_can_be_created() {
        let _world = World::new();
    }

    #[test]
    fn entity_lifecycle() {
        let mut world = World::new();
        let e = world.create_entity();
        assert!(world.alive_entities.contains(&e));
        world.delete_entity(e);
        assert!(!world.alive_entities.contains(&e));
    }

    #[derive(Debug, PartialEq)]
    struct TestComponent(i32);
    impl Component for TestComponent {}

    #[test]
    fn component_add_and_get() {
        let mut world = World::new();
        world.register_component::<TestComponent>();
        let e = world.create_entity();
        world.add_component(e, TestComponent(42));
        let c = world.get_component::<TestComponent>(e);
        assert_eq!(c.unwrap().0, 42);
    }
}
