#include "simulation/barnes_hut_simulation_with_collisions.h"

#include <algorithm>

#include "simulation/barnes_hut_simulation.h"
#include "simulation/naive_parallel_simulation.h"
#include <omp.h>

void BarnesHutSimulationWithCollisions::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void BarnesHutSimulationWithCollisions::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    // Tính toán lực và vị trí của các cơ thể (tương tự như trong simulate_epoch của BarnesHutSimulation)
    BarnesHutSimulation::simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);

    // Tìm và xử lý các va chạm
    find_collisions(universe);

    // Nếu yêu cầu vẽ hình, thực hiện vẽ
    if (create_intermediate_plots && (universe.current_simulation_epoch % plot_intermediate_epochs == 0)) {
        plotter.add_bodies_to_image(universe);
        plotter.write_and_clear();
    }
}

void BarnesHutSimulationWithCollisions::find_collisions(Universe& universe){
    const double COLLISION_DISTANCE_THRESHOLD = 100000000000.0; // 100,000,000 km

    // Duyệt qua tất cả các cặp thiên thể
    for (std::int32_t i = 0; i < universe.num_bodies; ++i) {
        for (std::int32_t j = i + 1; j < universe.num_bodies; ++j) {
            // Tính khoảng cách giữa cơ thể i và j
            Vector2d<double> direction = universe.positions[i] - universe.positions[j];
            double distance = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);

            // Kiểm tra xem khoảng cách có nhỏ hơn ngưỡng không
            if (distance < COLLISION_DISTANCE_THRESHOLD) {
                // Tính toán khối lượng và tốc độ
                double m1 = universe.weights[i];
                double m2 = universe.weights[j];
                Vector2d<double> v1 = universe.velocities[i];
                Vector2d<double> v2 = universe.velocities[j];

                // Chọn cơ thể nặng hơn và thực hiện va chạm
                if (m1 <= m2) {
                    // Cập nhật khối lượng và tốc độ của cơ thể thứ hai (nặng hơn)
                    universe.weights[j] = m1 + m2;
                    universe.velocities[j] = (v1.operator*(m1) + v2.operator*(m2)) / (m1 + m2);
                    // Cập nhật cơ thể thứ nhất (mất đi sau va chạm)
                    universe.weights[i] = 0;
                    universe.velocities[i] = Vector2d<double>(0, 0); // Tốc độ của cơ thể này sẽ trở thành 0
                } else {
                    // Cập nhật khối lượng và tốc độ của cơ thể thứ nhất (nặng hơn)
                    universe.weights[i] = m1 + m2;
                    universe.velocities[i] = (v1.operator*(m1) + v2.operator*(m2)) / (m1 + m2);
                    // Cập nhật cơ thể thứ hai (mất đi sau va chạm)
                    universe.weights[j] = 0;
                    universe.velocities[j] = Vector2d<double>(0, 0); // Tốc độ của cơ thể này sẽ trở thành 0
                }
            }
        }
    }

    for (std::int32_t k = universe.weights.size() - 1; k >= 0; --k) {
        if (universe.weights[k] == 0) {
            universe.weights.erase(universe.weights.begin() + k);
            universe.positions.erase(universe.positions.begin() + k);
            universe.velocities.erase(universe.velocities.begin() + k);
            universe.forces.erase(universe.forces.begin() + k);
            universe.num_bodies--;
        }
    }
}

void BarnesHutSimulationWithCollisions::find_collisions_parallel(Universe& universe){
    const double COLLISION_DISTANCE_THRESHOLD = 100000000000.0; // 100,000,000 km

    // Duyệt qua tất cả các cặp thiên thể
#pragma omp parallel shared(universe)
    {
        int number_of_threads = omp_get_num_threads();
        int work_for_thread = universe.num_bodies / number_of_threads +1;
        for (std::int32_t i = work_for_thread*omp_get_thread_num(); i < work_for_thread*(omp_get_thread_num() + 1); ++i) {
            if (i>=universe.num_bodies) {
                break;
            }
            for (std::int32_t j = i + 1; j < universe.num_bodies; ++j) {

                Vector2d<double> direction = universe.positions[i] - universe.positions[j];
                double distance = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);

                // Kiểm tra xem khoảng cách có nhỏ hơn ngưỡng không
                if (distance < COLLISION_DISTANCE_THRESHOLD) {
                    // Tính toán khối lượng và tốc độ
                    double m1 = universe.weights[i];
                    double m2 = universe.weights[j];
                    Vector2d<double> v1 = universe.velocities[i];
                    Vector2d<double> v2 = universe.velocities[j];

                    // Chọn cơ thể nặng hơn và thực hiện va chạm
                    if (m1 <= m2) {
                        // Cập nhật khối lượng và tốc độ của cơ thể thứ hai (nặng hơn)
                        universe.weights[j] = m1 + m2;
                        universe.velocities[j] = (v1.operator*(m1) + v2.operator*(m2)) / (m1 + m2);
                        // Cập nhật cơ thể thứ nhất (mất đi sau va chạm)
                        universe.weights[i] = 0;
                        universe.velocities[i] = Vector2d<double>(0, 0); // Tốc độ của cơ thể này sẽ trở thành 0
                    } else {
                        // Cập nhật khối lượng và tốc độ của cơ thể thứ nhất (nặng hơn)
                        universe.weights[i] = m1 + m2;
                        universe.velocities[i] = (v1.operator*(m1) + v2.operator*(m2)) / (m1 + m2);
                        // Cập nhật cơ thể thứ hai (mất đi sau va chạm)
                        universe.weights[j] = 0;
                        universe.velocities[j] = Vector2d<double>(0, 0); // Tốc độ của cơ thể này sẽ trở thành 0
                    }
                }
            }
        }
    }

    for (std::int32_t k = universe.weights.size() - 1; k >= 0; --k) {
        if (universe.weights[k] == 0) {
            universe.weights.erase(universe.weights.begin() + k);
            universe.positions.erase(universe.positions.begin() + k);
            universe.velocities.erase(universe.velocities.begin() + k);
            universe.forces.erase(universe.forces.begin() + k);
            universe.num_bodies--;
        }
    }
}
