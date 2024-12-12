//
// Created by QY.
//

#include "TSP.h"
#include <stdexcept>
#include <algorithm>
#include <stack>
#include <numeric>

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

    std::random_device rd;
    gen = std::mt19937(rd());

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

    tour = std::vector<int>(n, -1);
    best_tour = std::vector<int>(n, -1);
    in_best_tour = std::vector<std::vector<bool>>(n, std::vector<bool>(n, false));
    node_ptr = std::vector<Node*>(n, nullptr);
}

LKHSolver::~LKHSolver() {
    if (first_node == nullptr)
        return;

    Node* current = first_node;
    do {
        Node* next = current->suc;
        delete current;
        current = next;
    } while (current != first_node);

    first_node = nullptr;
}

void LKHSolver::solve() {
    create_candidate_set();
    if (!solved)
        solution = find_tour();
}

void LKHSolver::create_candidate_set() {
    int lower_bound = ascent();

    if (solved) {
        solution = lower_bound;
        generate_best_tour_from_one_tree();
        rearrange_tour();
        return;
    }

    minimum_one_tree(true);
    topological_sort();
    generate_candidates((int) excess * lower_bound);
}

int LKHSolver::ascent() {
    int best_one_tree_size, one_tree_size;
    unsigned int period = initial_period, p;
    bool initial_phase = true;

    one_tree_size = minimum_one_tree();
    best_one_tree_size = one_tree_size;

    // Check whether the minimal 1-tree is a valid tour
    if (all_zero(v)) {
        solved = true;
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
                solved = true;
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

int LKHSolver::minimum_one_tree(bool final_trial) {
    std::fill(v.begin(), v.end(), -2);

    int mst_size = minimum_spanning_tree(final_trial);
    find_second_closest_max_neighbor_leaf();

    ++v[edge_end_1];
    ++v[edge_end_2];
    if (final_trial) {
        in_best_tour[edge_end_1][edge_end_2] = true;
        in_best_tour[edge_end_2][edge_end_1] = true;
    }

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

int LKHSolver::minimum_spanning_tree(bool final_trial) {
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
            if (final_trial) {
                children[dad[next_vertex]].push_back(next_vertex);
                in_best_tour[next_vertex][dad[next_vertex]] = true;
                in_best_tour[dad[next_vertex]][next_vertex] = true;
            }
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

int LKHSolver::find_tour() {
    int best_cost = INT_MAX;
    int cost;

    for (int run = 0; run < run_times; ++run) {
        choose_initial_tour();
        generate_node_list();
        cost = lin_kernighan();
        if (cost < best_cost) {
            record_better_tour();
            best_cost = cost;
        }

        has_find_tour = true;
    }
    rearrange_tour();
    return best_cost;
}

void LKHSolver::choose_initial_tour() {
    std::vector<int> not_chosen_nodes(n);
    std::vector<bool> is_chosen(n, false);
    std::iota(not_chosen_nodes.begin(), not_chosen_nodes.end(), 0);
    std::uniform_int_distribution<> dis(0, n - 1);

    // Step 1: Choose a random starting node
    int i = dis(gen);
    is_chosen[i] = true;
    tour[0] = i;
    std::swap(not_chosen_nodes[i], not_chosen_nodes.back());
    not_chosen_nodes.pop_back();

    int tour_index = 1;
    int j, j_idx;

    // Step 2: Iteratively choose the next node
    while (!not_chosen_nodes.empty()) {
        std::vector<int> satisfied_nodes;

        // Try to find a node j that satisfies (a), (b), and (c)
        for (const int& candidate : candidate_nodes[i]) {
            if (candidate == -1) break;
            if (has_find_tour &&
                !is_chosen[candidate] &&
                candidate_alphas[i][candidate] == 0 &&
                in_best_tour[i][candidate]) {
                satisfied_nodes.push_back(candidate);
            }
        }

        // If no node satisfies all conditions, try to find a node satisfying (a)
        if (satisfied_nodes.empty()) {
            for (const int& candidate : candidate_nodes[i]) {
                if (candidate == -1) break;
                if (!is_chosen[candidate]) {
                    satisfied_nodes.push_back(candidate);
                }
            }
        }

        if (!satisfied_nodes.empty()) {
            dis = std::uniform_int_distribution<>(0, (int)satisfied_nodes.size() - 1);
            j = satisfied_nodes[dis(gen)];
            auto it = std::find(not_chosen_nodes.begin(), not_chosen_nodes.end(), j);
            j_idx = (int) std::distance(not_chosen_nodes.begin(), it);
        } else {
            dis = std::uniform_int_distribution<>(0, (int) not_chosen_nodes.size() - 1);
            j_idx = dis(gen);
            j = not_chosen_nodes[j_idx];
        }

        std::swap(not_chosen_nodes[j_idx], not_chosen_nodes.back());
        not_chosen_nodes.pop_back();

        tour[tour_index++] = j;
        is_chosen[j] = true;
        i = j;
    }
}

void LKHSolver::generate_node_list() {

    // Create nodes and link them to form the doubly linked list
    Node* prev_node = nullptr;
    for (int i = 0; i < n; ++i) {
        Node* curr_node = new Node{tour[i], i, prev_node, nullptr};
        node_ptr[tour[i]] = curr_node;
        if (prev_node != nullptr) {
            prev_node->suc = curr_node;
        } else {
            first_node = curr_node; // Set the first node
        }
        prev_node = curr_node;
    }

    // Complete the circular doubly linked list
    if (first_node && prev_node) {
        first_node->pred = prev_node;
        prev_node->suc = first_node;
    }
}

int LKHSolver::lin_kernighan() {
    Node *t1 = first_node, *t2;
    int cost = 0, gain;
    bool improved = true;

    do {
        cost += adj_node_cost(t1, t1->suc);
        t1 = t1->suc;
    } while (t1 != first_node);

    while (improved) {
        improved = false;
        do {
            t2 = t1->suc;
            gain = search_3opt(t1, t2);
            if (gain != -1) {
                store_tour();
                cost -= gain;
                improved = true;
            }
            t1 = t1->suc;
        } while (t1 != first_node);
    }

    return cost;
}

int LKHSolver::adj_node_cost(LKHSolver::Node *node1, LKHSolver::Node *node2) {
    return adj_matrix[node1->id][node2->id];
}

int LKHSolver::search_3opt(LKHSolver::Node *t1, LKHSolver::Node *t2) {
    Node *t3, *t4, *t5, *t6;
    int cost_left = 0, cost_right = adj_node_cost(t1, t2);

    for (const int& t2_candidate: candidate_nodes[t2->id]) {
        if (t2_candidate == -1) break;

        t3 = node_ptr[t2_candidate];
        if (t3 == t1 || t3 == t2->suc ||
                cost_left + adj_node_cost(t2, t3) >= cost_right) {
            continue;
        }
        t4 = t3->suc;
        cost_left += adj_node_cost(t2, t3);
        cost_right += adj_node_cost(t3, t4);

        for (const int& t4_candidate: candidate_nodes[t4->id]) {
            if (t4_candidate == -1) break;

            t5 = node_ptr[t4_candidate];
            if (!between(t5, t2, t3) || t5 == t3 ||
                    cost_left + adj_node_cost(t4, t5) >= cost_right) {
                continue;
            }
            t6 = t5->suc;
            cost_left += adj_node_cost(t4, t5);

            if (cost_left + adj_node_cost(t1, t6) < cost_right + adj_node_cost(t5, t6)) {
                make_3opt_move(t1, t2, t3, t4, t5, t6);
                return cost_right + adj_node_cost(t5, t6) - cost_left - adj_node_cost(t1, t6);
            }
            cost_left -= adj_node_cost(t4, t5);
        }
        cost_left -= adj_node_cost(t2, t3);
        cost_right -= adj_node_cost(t3, t4);
    }

    return -1;
}

void LKHSolver::record_better_tour() {
    best_tour = tour;
    std::fill(in_best_tour.begin(), in_best_tour.end(), std::vector<bool>(n, false));

    for (int i = 0; i < n - 1; ++i) {
        in_best_tour[tour[i]][tour[i + 1]] = true;
        in_best_tour[tour[i + 1]][tour[i]] = true;
    }
    in_best_tour[tour[0]][tour.back()] = true;
    in_best_tour[tour.back()][tour[0]] = true;
}

void LKHSolver::store_tour() {
    Node *node = first_node;

    do {
        tour[node->rank] = node->id;
    } while (node != first_node);
}

bool LKHSolver::between(LKHSolver::Node *t1, LKHSolver::Node *t2, LKHSolver::Node *t3) {
    // Use rank to determine the relative order in O(1) time.
    // Normalize ranks to handle circular comparisons.
    int rank1 = t1->rank;
    int rank2 = t2->rank;
    int rank3 = t3->rank;

    // Check if t1 is between t2 and t3 in the forward direction.
    // This works because in a circular list, (rank2 < rank1 < rank3)
    // or (rank3 < rank2 and rank1 wraps around) are the valid conditions.
    if ((rank2 <= rank1 && rank1 <= rank3) || (rank3 < rank2 && (rank1 >= rank2 || rank1 <= rank3))) {
        return true;
    }

    return false;
}

void LKHSolver::make_3opt_move(LKHSolver::Node *t1, LKHSolver::Node *t2, LKHSolver::Node *t3, LKHSolver::Node *t4,
                               LKHSolver::Node *t5, LKHSolver::Node *t6) {
    t1->suc = t6;
    t6->pred = t1;

    t3->suc = t2;
    t2->pred = t3;

    t5->suc = t4;
    t4->pred = t5;

    update_ranks();
}

void LKHSolver::update_ranks() {
    int r = 0;
    Node *node = first_node;

    do {
        node->rank = r++;
        node = node->suc;
    } while (node != first_node);
}

void LKHSolver::rearrange_tour() {
    // Check if the tour is empty
    if (best_tour.empty()) {
        return;
    }

    // Find the position of the start node in the tour
    auto it = std::find(best_tour.begin(), best_tour.end(), start);

    // Compute the index of the start node
    int start_index = (int) std::distance(best_tour.begin(), it);

    // Rearrange the tour to start from the specified node
    std::rotate(best_tour.begin(), best_tour.begin() + start_index, best_tour.end());
}

void LKHSolver::generate_best_tour_from_one_tree() {
    int right_count = 0, left_count = 0;
    int node = edge_end_1;

    while (node != -1) {
        best_tour[n - 1 - right_count] = node;
        node = dad[node];
        ++right_count;
    }

    node = edge_end_2;
    while (left_count + right_count < n) {
        best_tour[left_count] = node;
        node = dad[node];
        ++left_count;
    }
}
