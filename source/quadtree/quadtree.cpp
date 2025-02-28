#include "quadtree.h"

#include "quadtreeNode.h"
#include <set>
#include <algorithm>
#include <stdexcept>
#include <omp.h>
#include <numeric>

Quadtree::Quadtree(Universe& universe, BoundingBox bounding_box, std::int8_t construct_mode) {


    std::vector<std::int32_t> body_indices;
    for (std::int32_t i = 0; i < universe.num_bodies; ++i) {
        body_indices.push_back(i);
    }
    root = new QuadtreeNode(bounding_box);
    if (construct_mode == 0) {
        // Xây dựng tuần tự
        std::vector<QuadtreeNode*> nodes = construct(universe, bounding_box, body_indices);
        root->children = nodes;
    } else if (construct_mode == 1) {
        root->children = construct_task(universe, bounding_box, body_indices);
    } else if (construct_mode == 2) {
        root->children = construct_task_with_cutoff(universe, bounding_box, body_indices);;
    } else {
        throw std::invalid_argument("Invalid construct_mode");
    }
}

Quadtree::~Quadtree(){
        delete root;
}

void Quadtree::calculate_cumulative_masses(){
    root->calculate_node_cumulative_mass();
}

void Quadtree::calculate_center_of_mass(){
    root->calculate_node_center_of_mass();
}


std::vector<QuadtreeNode *> Quadtree::construct(Universe &universe, BoundingBox BB,
                                                std::vector<std::int32_t> body_indices) {
    if (body_indices.size() ==1) {
        QuadtreeNode* child_node = new QuadtreeNode(BB);
        child_node->body_identifier = body_indices[0];
        child_node->center_of_mass = universe.positions[body_indices[0]];
        child_node->cumulative_mass = universe.weights[body_indices[0]];
        child_node->center_of_mass_ready = true;
        child_node->cumulative_mass_ready = true;
        std::vector<QuadtreeNode*> child_nodeResult;
        child_nodeResult.push_back(child_node);
        return child_nodeResult;
    } else {
        std::vector<QuadtreeNode*> child_nodes;

        BoundingBox child_BBs[4] = {BB.get_quadrant(0),BB.get_quadrant(1),BB.get_quadrant(2),BB.get_quadrant(3)};
        for (BoundingBox& child_BB : child_BBs) {
            // Lọc các body nằm trong phần tư hiện tại
            std::vector<std::int32_t> bodies_in_child;
            for (std::int32_t body_index : body_indices) {
                Vector2d<double>& pos = universe.positions[body_index];
                if (child_BB.contains(pos)) {
                    bodies_in_child.push_back(body_index);
                }
            }
            if (bodies_in_child.empty()) {
                continue;
            }
            QuadtreeNode* child_node = new QuadtreeNode(child_BB);
            child_node->children = construct(universe, child_BB, bodies_in_child);
            child_nodes.push_back(child_node);
        }
        return child_nodes;
    }
}

std::vector<QuadtreeNode*> Quadtree::construct_task(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices) {
    if (body_indices.size() ==1) {
        QuadtreeNode* child_node = new QuadtreeNode(BB);
        child_node->body_identifier = body_indices[0];
        child_node->center_of_mass = universe.positions[body_indices[0]];
        child_node->cumulative_mass = universe.weights[body_indices[0]];
        child_node->center_of_mass_ready = true;
        child_node->cumulative_mass_ready = true;
        std::vector<QuadtreeNode*> child_nodeResult;
        child_nodeResult.push_back(child_node);
        return child_nodeResult;
    } else {
        std::vector<QuadtreeNode*> child_nodes;

        BoundingBox child_BBs[4] = {BB.get_quadrant(0),BB.get_quadrant(1),BB.get_quadrant(2),BB.get_quadrant(3)};
        #pragma omp parallel
        {
            #pragma omp single
            {
                for (BoundingBox& child_BB : child_BBs) {
                    #pragma omp task shared(child_nodes)
                    {
                        std::vector<std::int32_t> bodies_in_child;
                        for (std::int32_t body_index : body_indices) {
                            Vector2d<double>& pos = universe.positions[body_index];
                            if (child_BB.contains(pos)) {
                                bodies_in_child.push_back(body_index);
                            }
                        }
                        if (!bodies_in_child.empty()) {
                            QuadtreeNode* child_node = new QuadtreeNode(child_BB);
                            child_node->children = construct(universe, child_BB, bodies_in_child);
                            #pragma omp taskwait
                            child_nodes.push_back(child_node);
                        }
                    }
                }
            }
        }
        return child_nodes;
    }
}
std::vector<QuadtreeNode*> Quadtree::construct_task_with_cutoff(Universe& universe, BoundingBox& BB, std::vector<std::int32_t>& body_indices) {
     const int cutoff_threshold = 100;
     if (body_indices.size() == 1) {
        QuadtreeNode* child_node = new QuadtreeNode(BB);
        child_node->body_identifier = body_indices[0];
        child_node->center_of_mass = universe.positions[body_indices[0]];
        child_node->cumulative_mass = universe.weights[body_indices[0]];
        child_node->center_of_mass_ready = true;
        child_node->cumulative_mass_ready = true;
        std::vector<QuadtreeNode*> child_nodeResult;
        child_nodeResult.push_back(child_node);
        return child_nodeResult;
    }
        std::vector<QuadtreeNode*> child_nodes;
        BoundingBox child_BBs[4] = {BB.get_quadrant(0), BB.get_quadrant(1), BB.get_quadrant(2), BB.get_quadrant(3)};
        std::vector<std::vector<int32_t>> bodies_in_child(4);
        for (std::int32_t body_index : body_indices) {
            for (int i =0; i<4; i++) {
                Vector2d<double>& pos = universe.positions[body_index];
                if (child_BBs[i].contains(pos)) {
                    bodies_in_child[i].push_back(body_index);
                }
            }
        }
            #pragma omp parallel
            {
                #pragma omp single
                {
                    for (int i=0; i<4; i++) {
                    if (!bodies_in_child[i].empty()) {
                        #pragma omp task shared(child_nodes) final (bodies_in_child[i].size() < cutoff_threshold)
                        {
                            QuadtreeNode* child_node = new QuadtreeNode(child_BBs[i]);
                            if (bodies_in_child.size() <= cutoff_threshold) {
                                child_node->children = construct(universe, child_BBs[i], bodies_in_child[i]);
                            } else {
                                child_node->children = construct(universe, child_BBs[i], bodies_in_child[i]);
                            }
                            #pragma omp critical
                            child_nodes.push_back(child_node);
                        }
                    }
                    }
                }
            }
        return child_nodes;
}


std::vector<BoundingBox> Quadtree::get_bounding_boxes(QuadtreeNode* qtn){
    // traverse quadtree and collect bounding boxes
    std::vector<BoundingBox> result;
    // collect bounding boxes from children
    for(auto child: qtn->children){
        for(auto bb: get_bounding_boxes(child)){
            result.push_back(bb);
        }
    }
    result.push_back(qtn->bounding_box);
    return result;
}










