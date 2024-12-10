//
// Created by QY.
//

#ifndef TSP_H
#define TSP_H

#include <utility>
#include <vector>
#include <climits>

struct TspProblem {
    unsigned int n;
    int start;
    std::vector<std::vector<int>> adj_matrix;

    TspProblem(const std::vector<std::vector<int>>& adjacency_matrix, int s);
};

// A two-way linked list node
struct Node {
    int id;
    Node *pred, *suc;
};

class LKHSolver {
    unsigned int n;
    int start;
    std::vector<std::vector<int>> adj_matrix;
    int max_candidates;
    int run_times;

    Node *first_node;

public:
    explicit LKHSolver(TspProblem& p, int candidates = 5, int runs = 1);
    ~LKHSolver();

    void solve();

    void create_candidate_set();
    int ascent();
    [[nodiscard]] int minimum_one_tree() const;
    [[nodiscard]] int minimum_spanning_tree() const;

    int find_tour();
};

#endif