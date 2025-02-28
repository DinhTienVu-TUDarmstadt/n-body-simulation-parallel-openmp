#include "quadtreeNode.h"
#include "quadtree.h"

#include <iostream>


double QuadtreeNode::calculate_node_cumulative_mass(){
    if (children.empty()) {
        return cumulative_mass;
    } else {
        double total_mass = 0.0;

        for (auto& child : children) {
            total_mass += child->calculate_node_cumulative_mass();
        }

        return total_mass;
    }
}

QuadtreeNode::QuadtreeNode(BoundingBox arg_bounding_box)
    : bounding_box(arg_bounding_box), body_identifier(-1), cumulative_mass(0.0),
      center_of_mass_ready(false), cumulative_mass_ready(false) {
        children.clear();
}

QuadtreeNode::~QuadtreeNode(){
    for (auto child : children) {
        delete child;
    }
    children.clear();
}

Vector2d<double> QuadtreeNode::calculate_node_center_of_mass(){
    if (children.empty()) {
        return center_of_mass;
    } else {
        double total_mass = 0.0;
        Vector2d<double> weighted_position(0.0, 0.0);

        for (auto& child : children) {
            double child_mass = child->calculate_node_cumulative_mass();
            weighted_position = weighted_position + child->calculate_node_center_of_mass() * child_mass;
            total_mass += child_mass;
        }

        if (total_mass != 0) {
            center_of_mass = weighted_position / total_mass;
        } else {
            center_of_mass = Vector2d<double>(0.0, 0.0);
        }

        center_of_mass_ready = true;
        return center_of_mass;
    }
}