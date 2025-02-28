#include "simulation/naive_parallel_simulation.h"
#include "physics/gravitation.h"
#include "physics/mechanics.h"

#include <cmath>

void NaiveParallelSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void NaiveParallelSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    calculate_forces(universe);
    calculate_velocities(universe);
    calculate_positions(universe);
    universe.current_simulation_epoch++;
    if(create_intermediate_plots){
        if(universe.current_simulation_epoch % plot_intermediate_epochs == 0){
            plotter.add_bodies_to_image(universe);
            plotter.write_and_clear();
        }
    }
}


void NaiveParallelSimulation::calculate_forces(Universe& universe){
    std::size_t num_bodies = universe.num_bodies;
    universe.forces.clear();
    universe.forces.resize(num_bodies, Vector2d<double>(0, 0));

    // Song song hóa vòng lặp ngoài với OpenMP
#pragma omp parallel for
    for (size_t i = 0; i < num_bodies; i++) {
        Vector2d<double> force_sum(0, 0); // Lực tổng cho body thứ i

        for (size_t j = 0; j < num_bodies; j++) {
            if (i == j) continue; // Bỏ qua lực của chính nó

            Vector2d<double> direction = universe.positions[j] - universe.positions[i];
            double distance = sqrt(pow(direction[0], 2) + pow(direction[1], 2));

            if (distance > 0) {
                direction = direction / distance;
                double force_magnitude = gravitational_force(universe.weights[i], universe.weights[j], distance);
                Vector2d<double> force = direction * force_magnitude;
                force_sum = force_sum + force;
            }
        }

        // Cập nhật lực vào danh sách lực của universe
        universe.forces[i] = force_sum;
    }
}

void NaiveParallelSimulation::calculate_velocities(Universe& universe){
#pragma omp parallel for
    for (std::size_t i = 0; i < universe.num_bodies; i++) {
        // Lấy lực tác dụng trên body thứ i
        Vector2d<double> force = universe.forces[i];
        // Lấy khối lượng của body thứ i
        double mass = universe.weights[i];
        // Tính gia tốc
        Vector2d<double> acceleration = calculate_acceleration(force, mass);

        // Lấy vận tốc ban đầu và tính vận tốc mới
        Vector2d<double> initial_velocity = universe.velocities[i];
        Vector2d<double> new_velocity = calculate_velocity(initial_velocity, acceleration, 2.628e6);

        // Cập nhật vận tốc mới vào universe
        universe.velocities[i] = new_velocity;
    }
}

void NaiveParallelSimulation::calculate_positions(Universe& universe){
#pragma omp parallel for
    for (std::size_t i = 0; i < universe.num_bodies; ++i) {
        // Lấy vận tốc của body thứ i
        Vector2d<double> velocity = universe.velocities[i];

        // Tính toán khoảng cách di chuyển
        Vector2d<double> movement = velocity * 2.628e6;

        // Tính vị trí mới
        Vector2d<double> initial_position = universe.positions[i];
        Vector2d<double> new_position = initial_position + movement;

        // Cập nhật vị trí mới vào universe
        universe.positions[i] = new_position;
    }
}