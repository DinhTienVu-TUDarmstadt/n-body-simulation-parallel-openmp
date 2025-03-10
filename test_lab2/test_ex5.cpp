#include "test.h"

#include <exception>
#include <iostream>

#include "structures/universe.h"
#include "input_generator/input_generator.h"

#include "simulation/barnes_hut_simulation_with_collisions.h"

class Ex5Test : public LabTest {};

TEST_F(Ex5Test, test_five_a){
    Universe uni;

    // register collision pair
    uni.weights.push_back(100.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(200000000000.0, 200000000000.0));
    uni.velocities.push_back(Vector2d<double>(1000.0, 0.0));

    uni.weights.push_back(300.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(200000000000.0, 200500000000.0));
    uni.velocities.push_back(Vector2d<double>(600.0, -1000.0));

    uni.weights.push_back(600.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(200500000000.0, 200500000000.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    // register dummy bodies
    uni.weights.push_back(500.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(0.0, 0.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    uni.weights.push_back(700.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(-200500000000.0, 0.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    uni.weights.push_back(800.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(0, -200500000000.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    uni.weights.push_back(900.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(-200500000000.0, -200500000000.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    uni.num_bodies = 7;

    BarnesHutSimulationWithCollisions::find_collisions(uni);

    // check reduced amount of bodies
    ASSERT_EQ(uni.num_bodies, 5);
    ASSERT_EQ(uni.positions.size(), uni.num_bodies);
    ASSERT_EQ(uni.forces.size(), uni.num_bodies);
    ASSERT_EQ(uni.weights.size(), uni.num_bodies);
    ASSERT_EQ(uni.velocities.size(), uni.num_bodies);

    // check for collision of the correct bodies
    std::vector<double> valid_weights = {1000.0, 500.0, 700.0, 800.0, 900.0};

    std::int32_t valid_counter = 0;
    for(std::int32_t body_idx = 0; body_idx < uni.num_bodies; body_idx++){
        for(std::int32_t index = 0; index < valid_weights.size(); index++){
            if(valid_weights[index] == uni.weights[body_idx]){
                valid_counter++;
                valid_weights.erase(valid_weights.begin() + index);
                break;
            }
        }
    }  

    ASSERT_EQ(valid_counter, uni.num_bodies);

}


TEST_F(Ex5Test, test_five_b){
    Universe uni;

    // register collision pair
    uni.weights.push_back(100.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(200000000000.0, 200000000000.0));
    uni.velocities.push_back(Vector2d<double>(1000.0, 0.0));

    uni.weights.push_back(300.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(200000000000.0, 200500000000.0));
    uni.velocities.push_back(Vector2d<double>(600.0, -1000.0));

    uni.weights.push_back(600.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(200500000000.0, 200500000000.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    // register dummy bodies
    uni.weights.push_back(500.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(0.0, 0.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    uni.weights.push_back(700.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(-200500000000.0, 0.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    uni.weights.push_back(800.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(0, -200500000000.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    uni.weights.push_back(900.0);
    uni.forces.push_back(Vector2d<double>(0.0, 0.0));
    uni.positions.push_back(Vector2d<double>(-200500000000.0, -200500000000.0));
    uni.velocities.push_back(Vector2d<double>(6100.0, -1000.0));

    uni.num_bodies = 7;

    BarnesHutSimulationWithCollisions::find_collisions_parallel(uni);

    // check reduced amount of bodies
    ASSERT_EQ(uni.num_bodies, 5);
    ASSERT_EQ(uni.positions.size(), uni.num_bodies);
    ASSERT_EQ(uni.forces.size(), uni.num_bodies);
    ASSERT_EQ(uni.weights.size(), uni.num_bodies);
    ASSERT_EQ(uni.velocities.size(), uni.num_bodies);

    // check for collision of the correct bodies
    std::vector<double> valid_weights = {1000.0, 500.0, 700.0, 800.0, 900.0};

    std::int32_t valid_counter = 0;
    for(std::int32_t body_idx = 0; body_idx < uni.num_bodies; body_idx++){
        for(std::int32_t index = 0; index < valid_weights.size(); index++){
            if(valid_weights[index] == uni.weights[body_idx]){
                valid_counter++;
                valid_weights.erase(valid_weights.begin() + index);
                break;
            }
        }
    }

    ASSERT_EQ(valid_counter, uni.num_bodies);

}

TEST_F(Ex5Test, test_five_b_runtime){
    // initialize
    Universe uni;
    InputGenerator::create_random_universe(100, uni);
    BoundingBox BB = uni.get_bounding_box();
    // construct quadtree sequentially
    // validate quadtree
    auto start1 = std::chrono::high_resolution_clock::now();

    // Code block to measure
    BarnesHutSimulationWithCollisions::find_collisions(uni);

    // End time
    auto end1 = std::chrono::high_resolution_clock::now();

    // Calculate runtime in nanoseconds
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);

    auto start2 = std::chrono::high_resolution_clock::now();

    // Code block to measure
    BarnesHutSimulationWithCollisions::find_collisions_parallel(uni);

    // End time
    auto end2 = std::chrono::high_resolution_clock::now();

    // Calculate runtime in nanoseconds
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2- start2);


    std::cout << "Runtime: " << duration1.count() << " microseconds" << std::endl;
    std::cout << "Runtime: " << duration2.count() << " microseconds" << std::endl;
    ASSERT_TRUE(duration2<=duration1);
}




