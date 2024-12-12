//
// Created by QY.
//

#ifndef TSP_H
#define TSP_H

#include <utility>
#include <vector>
#include <climits>
#include <iostream>
#include <random>

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

    std::mt19937 gen;

    // Default parameters
    int initial_step_size = 1;
    unsigned int initial_period;  // Default by max(n / 2, 100)
    float excess = 1.0;

    // Vectors that will be used for searching candidates
    std::vector<int> dad;  // store parents of the nodes
    std::vector<std::vector<int>> children;  // store children of the nodes
    std::vector<int> v;  // store degrees of the nodes, v[i] = degree[i] - 2
    std::vector<bool> on_tree;  // Track nodes included in MST
    std::vector<int> topological_order;  // Track the topological order of the nodes
    std::vector<int> key;  // key values to select minimum weight edge
    std::vector<int> pi;  // pi value for each node
    std::vector<int> best_pi;  // best pi value for each node
    std::vector<std::vector<int>> candidate_nodes;
    std::vector<std::vector<int>> candidate_alphas;

    // The additional edge from MST to 1-tree
    int edge_end_1 = 0;
    int edge_end_2 = 0;

    struct Node {
        int id;  // the identification number of the node
        int rank;  // Rank gives the ordinal number of the node in the tour
        Node *pred, *suc;  // point to the predecessor node and the successor node of the tour, respectively
    };

    // Vectors that will be used for finding best tour
    std::vector<int> tour;
    std::vector<std::vector<bool>> in_best_tour;
    std::vector<Node*> node_ptr;
    bool has_find_tour = false;

    // Variables for finding best tour

    Node* first_node = nullptr;
    bool solved = false;

    // Part 1: Create candidate set
    void create_candidate_set();
    int ascent();
    int minimum_one_tree(bool final_trial = false);
    [[nodiscard]] int minimum_spanning_tree(bool final_trial = false);
    void generate_candidates(int max_alpha);
    void insert_candidate(int from, int to, int alpha);

    static bool all_zero(std::vector<int>& vec);
    int best_pi_edge(int& from, int& to) const;
    void find_second_closest_max_neighbor_leaf();
    void topological_sort();

    // Part 2: Solve the problem
    int find_tour();
    void choose_initial_tour();
    void generate_node_list();

    int lin_kernighan();
    int adj_node_cost(Node *node1, Node *node2);
    int search_3opt(Node *t1, Node *t2);
    void store_tour();
    static bool between(Node *t1, Node *t2, Node *t3);
    void make_3opt_move(Node *t1, Node *t2, Node *t3, Node *t4, Node *t5, Node *t6);
    void update_ranks();

    void record_better_tour();

public:
    explicit LKHSolver(TspProblem& p, int candidates = 5, int runs = 1);
    ~LKHSolver();

    void solve();
    int solution = -1;
    std::vector<int> best_tour;
};

#endif