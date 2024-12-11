//
// Created by QY.
//

#include "TSP.h"
#include <stdexcept>
#include <algorithm>
#include <stack>

TspProblem::TspProblem(const std::vector<std::vector<int>>& adjacency_matrix, int s)
        : n((int) adjacency_matrix.size()), start(s), adj_matrix(adjacency_matrix) {
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

    candidate_nodes = std::vector<std::vector<int>>(n, std::vector<int>(max_candidates, -1));
    candidate_alphas = std::vector<std::vector<int>>(n, std::vector<int>(max_candidates, INT_MAX));

    initial_period = std::max(n / 2, 100);

    dad = std::vector<int>(n, -1);
    children = std::vector<std::vector<int>>(n, std::vector<int>());
    v = std::vector<int>(n, 0);
    on_tree = std::vector<bool>(n, false);
    topological_order = std::vector<int>(n, 0);
    key = std::vector<int>(n, INT_MAX);
    pi = std::vector<int>(n, 0);
    best_pi = std::vector<int>(n, 0);
}

LKHSolver::~LKHSolver() = default;

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

    minimum_one_tree(true);
    topological_sort();
    generate_candidates((int) excess * lower_bound);
}

int LKHSolver::find_tour() {
    return 0;
}

int LKHSolver::ascent() {
    int best_one_tree_size, one_tree_size;
    unsigned int period = initial_period, p;
    bool initial_phase = true;

    one_tree_size = minimum_one_tree();
    best_one_tree_size = one_tree_size;

    // Check whether the minimal 1-tree is a valid tour

    if (all_zero(v)) {
        return one_tree_size;
    }

    std::vector<int> last_v = v;

    // Perform sub-gradient descent
    for (int t = initial_step_size; t > 0; period /= 2, t /= 2) {
        for (p = 1; t > 0 && p <= period; ++p) {
            for (int i = 0; i < n; ++i) {
                if (v[i] != 0)
                    pi[i] += t * (7 * v[i] + 3 * last_v[i]) / 10;
                last_v[i] = v[i];
            }

            one_tree_size = minimum_one_tree();
            if (all_zero(v)) {
                return one_tree_size;
            }

            if (one_tree_size > best_one_tree_size) {
                best_one_tree_size = one_tree_size;
                for (int i = 0; i < n; ++i)
                    best_pi[i] = pi[i];
                if (initial_phase)
                    t *= 2;
                if (p == period)
                    period *= 2;
            } else if (initial_phase && p > initial_period / 2) {
                initial_phase = false;
                p = 0;
                t = 3 * t / 4;
            }
        }
    }

    return best_one_tree_size;
}

int LKHSolver::minimum_one_tree(bool maintain_children) {
    std::fill(v.begin(), v.end(), -2);

    int mst_size = minimum_spanning_tree(maintain_children);
    find_second_closest_max_neighbor_leaf();

    ++v[edge_end_1];
    ++v[edge_end_2];

    int sum_pi = 0;
    for (const int& pi_value: pi) {
        sum_pi += pi_value;
    }

    return mst_size + adj_matrix[edge_end_1][edge_end_2] + pi[edge_end_1] + pi[edge_end_2] - 2 * sum_pi;
}

void LKHSolver::find_second_closest_max_neighbor_leaf() {
    int max_second_closest_distance = INT_MIN;

    for (int i = 0; i < n; ++i) {
        if (v[i] != -1) continue;

        int first_min = INT_MAX;
        int second_min = INT_MAX;
        int first_min_id = -1;
        int second_min_id = -1;

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            int distance = adj_matrix[i][j] + pi[i] + pi[j];
            if (distance < first_min) {
                second_min = first_min;
                second_min_id = first_min_id;
                first_min = distance;
                first_min_id = j;
            } else if (distance < second_min) {
                second_min = distance;
                second_min_id = j;
            }
        }

        // Update the max second-closest distance and the node id
        if (second_min > max_second_closest_distance) {
            max_second_closest_distance = second_min;
            edge_end_1 = i;
            edge_end_2 = second_min_id;
        }
    }
}

int LKHSolver::minimum_spanning_tree(bool maintain_children) {
    // TODO: May use a heap to speed up when the graph is sparse
    std::fill(key.begin(), key.end(), INT_MAX);
    std::fill(on_tree.begin(), on_tree.end(), false);
    dad[0] = -1;
    key[0] = 0;

    int mst_size = 0;
    for (int count = 0; count < n; ++count) {

        // Find the vertex with the minimum key value, not yet included in MST
        int next_vertex = -1;
        int min_key = INT_MAX;

        for (int i = 0; i < n; ++i) {
            if (!on_tree[i] && key[i] < min_key) {
                next_vertex = i;
                min_key = key[i];
            }
        }

        // Include this vertex in the MST
        on_tree[next_vertex] = true;
        mst_size += min_key;

        if (count != 0) {
            ++v[next_vertex];
            ++v[dad[next_vertex]];
            if (maintain_children)
                children[dad[next_vertex]].push_back(next_vertex);
        }

        // Update the key values of adjacent vertices
        for (int neighbor = 0; neighbor < n; ++neighbor) {
            int dist = adj_matrix[next_vertex][neighbor] + pi[next_vertex] + pi[neighbor];
            if (neighbor != next_vertex && !on_tree[neighbor] && dist < key[neighbor]) {
                key[neighbor] = dist;
                dad[neighbor] = next_vertex;
            }
        }
    }

    return mst_size;
}

void LKHSolver::topological_sort() {
    std::stack<int> dfs_stack;
    std::vector<bool> visited(n, false);
    int curr_node = edge_end_1;
    int last_node = -1;
    int tmp_dad;
    int count = 0;

    // Perform DFS
    while (curr_node != -1) {
        for (const int& child : children[curr_node]) {
            if (!visited[child]) {
                dad[child] = curr_node;
                dfs_stack.push(child);
            }
        }
        // Update topological order
        topological_order[count] = curr_node;
        ++count;
        visited[curr_node] = true;

        // Update the new parent relationship, O(1) space complexity
        tmp_dad = dad[curr_node];
        dad[curr_node] = last_node;
        last_node = curr_node;
        curr_node = tmp_dad;

        while (!dfs_stack.empty()) {
            int node = dfs_stack.top();
            dfs_stack.pop();
            topological_order[count] = node;
            ++count;

            // Push children of the current node to stack
            for (int child : children[node]) {
                dad[child] = node;
                dfs_stack.push(child);
            }
        }
    }
}

void LKHSolver::generate_candidates(int max_alpha) {
    int alpha;
    int father;
    std::vector<int> mark = std::vector<int>(n, 0);
    std::vector<int> beta = std::vector<int>(n, 0);
    int edge_one_tree = best_pi_edge(edge_end_1, edge_end_2);

    for (int f = 0; f < n; ++f) {
        int from = topological_order[f];
        if (from != edge_end_1) {
            beta[from] = INT_MIN;
            for (int to = from; dad[to] != -1; to = dad[to]) {
                father = dad[to];
                beta[father] = std::max(beta[to], best_pi_edge(father, to));
                mark[father] = from;
            }
        }

        for (int t = 0; t < n; ++t) {
            if (f == t)
                continue;
            int to = topological_order[t];

            if (from == edge_end_1) {
                alpha = to == edge_end_2 || dad[to] == edge_end_1 ? 0 :
                        best_pi_edge(from, to) - edge_one_tree;
            } else if (to == edge_end_1) {
                alpha = from == edge_end_2 || dad[from] == edge_end_1 ? 0 :
                        best_pi_edge(from, to) - edge_one_tree;
            } else {
                if (mark[to] != from) {
                    father = dad[to];
                    beta[to] = std::max(beta[father], best_pi_edge(father, to));
                }
                alpha = best_pi_edge(from, to) - beta[to];
            }

            if (alpha <= max_alpha)
                insert_candidate(from, to, alpha);
        }
    }
}

void LKHSolver::insert_candidate(int from, int to, int alpha) {
    for (int i = 0; i < max_candidates; ++i) {
        if (candidate_nodes[from][i] == -1) {
            candidate_nodes[from][i] = to;
            candidate_alphas[from][i] = alpha;
            return;
        }

        if (alpha < candidate_alphas[from][i]) {
            for (int j = 4; j >= i + 1; --j) {
                candidate_alphas[from][j] = candidate_alphas[from][j - 1];
                candidate_nodes[from][j] = candidate_nodes[from][j - 1];
            }
            candidate_nodes[from][i] = to;
            candidate_alphas[from][i] = alpha;
            return;
        }
    }
}


bool LKHSolver::all_zero(std::vector<int> &vec) {
    return std::all_of(vec.begin(), vec.end(), [](int v_value) { return v_value == 0; });
}

int LKHSolver::best_pi_edge(int &from, int &to) const {
    return adj_matrix[from][to] + best_pi[from] + best_pi[to];
}

