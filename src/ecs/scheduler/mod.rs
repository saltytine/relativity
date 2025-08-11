//! ecs scheduler and dispatcher
//!
//! runs systems in order, supports event-driven execution

use super::{System, World};


/// event queue for ecs events
pub struct EventQueue {
    events: Vec<Box<dyn std::any::Any + Send + Sync>>,
}

impl EventQueue {
    pub fn new() -> Self {
        Self { events: Vec::new() }
    }
    pub fn push<E: 'static + Send + Sync>(&mut self, event: E) {
        self.events.push(Box::new(event));
    }
    pub fn drain(&mut self) -> Vec<Box<dyn std::any::Any + Send + Sync>> {
        std::mem::take(&mut self.events)
    }
}

/// scheduler: runs systems in order, supports event dispatch
pub struct Scheduler {
    systems: Vec<Box<dyn System>>,
}

impl Scheduler {
    pub fn new() -> Self {
        Self { systems: Vec::new() }
    }
    pub fn add_system<S: System + 'static>(&mut self, system: S) {
        self.systems.push(Box::new(system));
    }
    pub fn run(&mut self, world: &mut World) {
        for sys in self.systems.iter_mut() {
            sys.run(world);
        }
    }
}

/// dispatcher: runs systems, can process events (future extension)
pub struct Dispatcher {
    scheduler: Scheduler,
    event_queue: EventQueue,
}

impl Dispatcher {
    pub fn new() -> Self {
        Self {
            scheduler: Scheduler::new(),
            event_queue: EventQueue::new(),
        }
    }
    pub fn add_system<S: System + 'static>(&mut self, system: S) {
        self.scheduler.add_system(system);
    }
    pub fn push_event<E: 'static + Send + Sync>(&mut self, event: E) {
        self.event_queue.push(event);
    }
    pub fn run(&mut self, world: &mut World) {
        self.scheduler.run(world);
        // event processing can be added here
    }
}
