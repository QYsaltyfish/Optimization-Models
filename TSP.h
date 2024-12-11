//
// Created by QY.
//

#ifndef TSP_H
#define TSP_H

#include <utility>
#include <vector>
#include <climits>
#include <iostream>

struct TspProblem {
    int n;
    int start;
    std::vector<std::vector<int>> adj_matrix;

    TspProblem(const std::vector<std::vector<int>>& adjacency_matrix, int s);
};

class LKHSolver {
    int n;
    int start;
    std::vector<std::vector<int>> adj_matrix;
    int max_candidates;
    int run_times;

    std::vector<std::vector<int>> candidate_nodes;
    std::vector<std::vector<int>> candidate_alphas;

    // Default parameters
    int initial_step_size = 1;
    unsigned int initial_period;  // Default by max(n / 2, 100)
    float excess = 1.0;

    // Vectors that will be frequently used for searching candidates
    std::vector<int> dad;  // store parents of the nodes
    std::vector<std::vector<int>> children;  // store children of the nodes
    std::vector<int> v;  // store degrees of the nodes, v[i] = degree[i] - 2
    std::vector<bool> on_tree;  // Track nodes included in MST
    std::vector<int> topological_order;  // Track the topological order of the nodes
    std::vector<int> key;  // key values to select minimum weight edge
    std::vector<int> pi;  // pi value for each node
    std::vector<int> best_pi;  // best pi value for each node

    // The additional edge from MST to 1-tree
    int edge_end_1 = 0;
    int edge_end_2 = 0;

    bool solved = false;
    int solution = -1;

    static bool all_zero(std::vector<int>& vec);
    int best_pi_edge(int& from, int& to) const;
    void find_second_closest_max_neighbor_leaf();
    void topological_sort();
public:
    explicit LKHSolver(TspProblem& p, int candidates = 5, int runs = 1);
    ~LKHSolver();

    void solve();

    void create_candidate_set();
    int ascent();
    int minimum_one_tree(bool maintain_children = false);
    [[nodiscard]] int minimum_spanning_tree(bool maintain_children = false);
    void generate_candidates(int max_alpha);
    void insert_candidate(int from, int to, int alpha);

    int find_tour();
};

#endif