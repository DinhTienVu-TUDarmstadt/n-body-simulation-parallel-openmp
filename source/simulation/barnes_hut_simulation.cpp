#include "simulation/barnes_hut_simulation.h"
#include "simulation/naive_parallel_simulation.h"
#include "physics/gravitation.h"
#include "physics/mechanics.h"
#include "plotting/plotter.h"

#include <cmath>
#include <functional>

void BarnesHutSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void BarnesHutSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    BoundingBox universe_bb = universe.get_bounding_box();
    std::vector<std::int32_t> all_body_indices;
    for (std::int32_t i = 0; i < universe.num_bodies; i++) {
        all_body_indices.push_back(i);
    }

    Quadtree quadtree(universe, universe_bb, 0);  // Quadtree với mức độ 0
    quadtree.root->children = quadtree.construct(universe, universe_bb, all_body_indices);  // Xây dựng Quadtree

    std::vector<QuadtreeNode*> queue = {quadtree.root};
    while (!queue.empty()) {
        auto current_node = queue.back();
        queue.pop_back();

        current_node->calculate_node_cumulative_mass();
        current_node->calculate_node_center_of_mass();

        for (auto& child : current_node->children) {
            queue.push_back(child);
        }
    }

    calculate_forces(universe,quadtree);

    NaiveParallelSimulation::calculate_velocities(universe);
    NaiveParallelSimulation::calculate_positions(universe);

    universe.current_simulation_epoch++;

    if (create_intermediate_plots && (universe.current_simulation_epoch % plot_intermediate_epochs == 0)) {
        plotter.add_bodies_to_image(universe);
        plotter.write_and_clear();
    }
}



void BarnesHutSimulation::get_relevant_nodes(Universe& universe, Quadtree& quadtree, std::vector<QuadtreeNode*>& relevant_nodes, Vector2d<double>& body_position, std::int32_t body_index, double threshold_theta) {
    std::vector<QuadtreeNode*> all_vectors;
    all_vectors.push_back(quadtree.root);
    while (!all_vectors.empty()) {
        auto current_node = all_vectors.back();
        all_vectors.pop_back();

        double diagonal = current_node->bounding_box.get_diagonal();
        Vector2d<double> direction = current_node->center_of_mass - body_position;
        double distance = sqrt(pow(direction[0], 2) + pow(direction[1], 2));
        double theta = diagonal / distance;
        if (current_node->body_identifier == body_index) {
            continue;
        }
        if (current_node->bounding_box.contains(body_position)) {
            for (QuadtreeNode* child : current_node->children) {
                all_vectors.push_back(child);
            }
        }
        else if (theta <= threshold_theta) {
            relevant_nodes.push_back(current_node);
        }
        else if (current_node->body_identifier != -1) {
            relevant_nodes.push_back(current_node);
        } else {
            for (QuadtreeNode* child : current_node->children) {
                all_vectors.push_back(child);
            }
        }
    }
}

void BarnesHutSimulation::calculate_forces(Universe& universe, Quadtree& quadtree){
    const double threshold_theta = 0.2;
    const double G = 6.67430e-11;

#pragma omp parallel for
    for (std::int32_t body_index = 0; body_index < universe.num_bodies; body_index++) {
        // Vị trí và các thông tin của cơ thể
        Vector2d<double>& body_position = universe.positions[body_index];
        Vector2d<double>& body_velocity = universe.velocities[body_index];
        double body_force_x = universe.forces[body_index][0];
        double body_force_y = universe.forces[body_index][1];

        body_force_x = 0;
        body_force_y = 0;

        std::vector<QuadtreeNode*> relevant_nodes;
        get_relevant_nodes(universe, quadtree, relevant_nodes, body_position, body_index, threshold_theta);

        for (QuadtreeNode* node : relevant_nodes) {
            if (node->center_of_mass_ready) {
                Vector2d<double> direction = node->center_of_mass - body_position;
                double distance = sqrt(pow(direction[0], 2) + pow(direction[1], 2));

                if (distance > 0) {
                    double force_magnitude = (G * universe.weights[body_index] * node->cumulative_mass) / (distance * distance);
                    direction = direction / distance;

                    body_force_x += force_magnitude * direction[0];
                    body_force_y += force_magnitude * direction[1];
                }
            }
        }
    }
}