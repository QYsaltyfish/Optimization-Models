//
// Created by QY.
//

#include "TSP.h"
#include <stdexcept>

TspProblem::TspProblem(const std::vector<std::vector<int>>& adjacency_matrix, int s)
        : n(adjacency_matrix.size()), start(s), adj_matrix(adjacency_matrix) {
    if (adjacency_matrix.size() != adjacency_matrix[0].size()) {
        throw std::invalid_argument("The adjacency matrix is not square.");
    }
}


LKHSolver::LKHSolver(TspProblem &p, int candidates, int runs) {
    n = p.n;
    start = p.start;
    adj_matrix = p.adj_matrix;

    max_candidates = candidates;
    run_times = runs;

    // Initialize the linked list
    Node *previous_node = nullptr;
    for (int i = 0; i < n; ++i) {
        Node *new_node = new Node();
        new_node->id = i;
        new_node->pred = previous_node;
        new_node->suc = nullptr;

        if (i == 0) {
            first_node = new_node;
        }

        if (previous_node != nullptr) {
            previous_node->suc = new_node;
        }
        previous_node = new_node;
    }

    if (previous_node != nullptr && first_node != nullptr) {
        previous_node->suc = first_node;
        first_node->pred = previous_node;
    }
}

LKHSolver::~LKHSolver() {
    Node *curr_node = first_node;
    Node *next_node = curr_node->suc;
    while (next_node != first_node) {
        delete curr_node;
        curr_node = next_node;
        next_node = next_node->suc;
    }
}

void LKHSolver::solve() {
    create_candidate_set();
    int best_cost = INT_MAX;

    for (int run = 0; run < run_times; ++run) {
        int cost = find_tour();
        best_cost = cost < best_cost ? cost : best_cost;
    }
}

void LKHSolver::create_candidate_set() {
    int lower_bound = ascent();
}

int LKHSolver::find_tour() {
    return 0;
}

int LKHSolver::ascent() {
    int best_w, w;
    int period = 1, p = 0;
    bool initial_phase = true;

    return 0;
}

int LKHSolver::minimum_one_tree() const {
    return 0;
}

int LKHSolver::minimum_spanning_tree() const {
    return 0;
}

