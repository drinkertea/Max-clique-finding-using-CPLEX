#include "solution.h"

#include <algorithm>
#include <random>
#include <iostream>

using namespace std;

mt19937 generator; // Mersenne Twister 19937 generator

SolutionMIS::SolutionMIS(const Graph* g, const std::vector<int>& weights) :
    g_(g),
    solution_(g_->GetSize()),
    solution_size_(0),
    free_size_(g_->GetSize()),
    tightness_(g_->GetSize(), 0),
    position_(g_->GetSize()),
    mu_(g_->GetSize()),
    weight_(0),
    weights(&weights)
{
    for (int idx = 0; idx < g_->GetSize(); idx++) {
        position_[idx] = idx;
        solution_[idx] = idx;
        mu_[idx] = weights[idx];
    }
    generator.seed(42);
} // SolutionMIS::SolutionMIS(const Graph *g)

void SolutionMIS::moveFreeToSolutionPartition(const int v)
{
    assert(v < g_->GetSize());

    // current position of v in the solution_ vector
    int pos_v = position_[v];

    // new position of v in the solution_ vector
    int new_pos_v = solution_size_;

    // first vertex of the second partition
    int j = solution_[solution_size_];

    // ensures v is in the free partition of the solution vector
    assert((solution_size_ <= pos_v) && (solution_size_ + free_size_ > pos_v));

    // swap v with the first vertex of the second partition
    swap(solution_[pos_v], solution_[new_pos_v]);
    position_[v] = new_pos_v;
    position_[j] = pos_v;

    // change the boundary between the blocks to make v the last vertex of the
    // first partition
    solution_size_++;
    free_size_--;
} // void SolutionMIS::moveFreeToSolutionPartition(const int v)

void SolutionMIS::moveFreeToNonFreePartition(const int v)
{
    assert(v < g_->GetSize());

    // current position of v in the solution vector
    int pos_v = position_[v];

    // new position of v in the solution vector
    int new_pos_v = solution_size_ + free_size_ - 1;

    // last vertex of the second partition
    int j = solution_[solution_size_ + free_size_ - 1];

    // ensures v is in the free partition of the solution vector
    assert((solution_size_ <= pos_v) && (solution_size_ + free_size_ > pos_v));

    // swap v with the last vertex of the second partition
    swap(solution_[pos_v], solution_[new_pos_v]);
    position_[v] = new_pos_v;
    position_[j] = pos_v;

    // change the boundary between the blocks to make v the last vertex of the
    // second partition
    free_size_--;
} // void SolutionMIS::moveFreeToNonFreePartition(const int v)

void SolutionMIS::moveSolutionToFreePartition(const int v)
{
    assert(v < g_->GetSize());

    // current position of v in the solution vector
    int pos_v = position_[v];

    // new position of v in the solution vector
    int new_pos_v = solution_size_ - 1;

    // last vertex of the first partition
    int j = solution_[solution_size_ - 1];

    // ensures v is in the solution partition of the solution vector
    assert(pos_v < solution_size_);

    // swap v with the last vertex of the second partition
    swap(solution_[pos_v], solution_[new_pos_v]);
    position_[v] = new_pos_v;
    position_[j] = pos_v;

    // change the boundary between the blocks to make v the first vertex of the
    // second partition
    solution_size_--;
    free_size_++;
} // void SolutionMIS::moveSolutionToFreePartition(const int v)

void SolutionMIS::moveNonFreeToFreePartition(const int v)
{
    assert(v < g_->GetSize());

    // current position of v in the solution vector
    int pos_v = position_[v];

    // new position of v in the solution vector
    int new_pos_v = solution_size_ + free_size_;

    // first vertex of the third partition
    int j = solution_[solution_size_ + free_size_];

    // ensures v is in the non free partition of the solution vector
    assert(pos_v >= solution_size_ + free_size_);

    // swap v with the last vertex of the second partition
    swap(solution_[pos_v], solution_[new_pos_v]);
    position_[v] = new_pos_v;
    position_[j] = pos_v;

    // change the boundary between the blocks to make v the last vertex of the
    // second partition
    free_size_++;
} // void SolutionMIS::moveNonFreeToFreePartition(const int v)

void SolutionMIS::addVertex(const int v)
{
    int weight_v = (*weights)[v];
    weight_ += weight_v;

    moveFreeToSolutionPartition(v);

    const vector<int>& adj_l = g_->GetAdj(v);

    for (int neighbor : adj_l) {
        // increase the tighness of each neighbor by one
        tightness_[neighbor]++;

        mu_[neighbor] -= weight_v;

        // if the neighbor is in the free partition, move to non free partition
        int neighbor_pos = position_[neighbor];
        if ((solution_size_ <= neighbor_pos) && (solution_size_ + free_size_ > neighbor_pos)) {
            moveFreeToNonFreePartition(neighbor);
        }
    }
} // void SolutionMIS::addVertex(const int v)

void SolutionMIS::removeVertex(const int v)
{
    int weight_v = (*weights)[v];
    weight_ -= weight_v;

    moveSolutionToFreePartition(v);

    const vector<int>& adj_l = g_->GetAdj(v);

    for (int neighbor : adj_l) {
        tightness_[neighbor]--;

        mu_[neighbor] += weight_v;

        // if the neighbor becomes free
        if (tightness_[neighbor] == 0) {
            moveNonFreeToFreePartition(neighbor);
        }
    }
} // void SolutionMIS::removeVertex(const int v)

bool SolutionMIS::integrityCheck() const
{
    for (int idx = 0; idx < solution_size_; idx++) {
        int vertex = solution_[idx];

        if (tightness_[vertex] > 0) {
            return false;
        }

        for (int neighbor : g_->GetAdj(vertex)) {
            if (find(solution_.begin(), solution_.begin() + solution_size_, neighbor)
                != solution_.begin() + solution_size_) {
                return false;
            }
        }
    }

    for (int idx = solution_size_; idx < solution_size_ + free_size_; idx++) {
        int vertex = solution_[idx];
        if (tightness_[vertex] > 0) {
            return false;
        }
    }

    for (int idx = solution_size_ + free_size_; idx < g_->GetSize(); idx++) {
        int vertex = solution_[idx];
        if (tightness_[vertex] == 0) {
            return false;
        }
    }

    return true;
} // bool SolutionMIS::integrityCheck() const

void SolutionMIS::addRandomVertex()
{
    assert(!isMaximal());

    // generate a random number between [0, free_size_ - 1]
    uniform_int_distribution<int> distribution(0, free_size_ - 1);
    int free_pos = distribution(generator);

    int vertex = solution_[solution_size_ + free_pos];

    addVertex(vertex);
} // void SolutionMIS::addRandomVertex()

bool SolutionMIS::omegaImprovement()
{
    for (int idx = g_->GetSize() - 1; idx >= solution_size_; idx--) {
        int v = solution_[idx];
        if (mu_[v] > 0) {
            for (int neighbor : g_->GetAdj(v)) {
                if (position_[neighbor] < solution_size_) {
                    removeVertex(neighbor);
                }
            }
            addVertex(v);
            return true;
        }
    }

    return false;
} // bool SolutionMIS::swapImprovement()

bool SolutionMIS::twoImprovement()
{
    assert(isMaximal());

    for (int idx = 0; idx < solution_size_; idx++) {
        // the candidate for removal
        int x = solution_[idx];

        // sorted list of 1-tight nighbors of x
        vector<int> onetight_list;

        // build the list of 1-tight nighbors of x
        for (int neighbor : g_->GetAdj(x)) {
            if (tightness_[neighbor] == 1) {
                onetight_list.push_back(neighbor);
            }
        }
        assert(is_sorted(onetight_list.begin(), onetight_list.end()));

        // if x has fewer than two 1-tight neighors we are done with x
        if (onetight_list.size() < 2) continue;

        int x_weight = (*weights)[x];

        // attempt to find in onetight_list a pair {v, w} such that there
        // is no edge between v and w
        for (int v : onetight_list) {

            // stores the sorted list of v nighbors
            vector<int> v_neighbors(g_->GetAdj(v));
            assert(is_sorted(v_neighbors.begin(), v_neighbors.end()));

            // check if there is a vertex w in onetight_list (besides v) that
            // does not belong to v_neighbors. since both onetight_list and v_neighbors
            // are sorted, this can be done by traversing both lists in tandem.
            size_t i_idx = 0, j_idx = 0;
            while (i_idx < v_neighbors.size()
                && j_idx < onetight_list.size()) {
                if (onetight_list[j_idx] == v) {
                    j_idx++;
                    continue;
                }
                else if (v_neighbors[i_idx] < onetight_list[j_idx]) {
                    i_idx++;
                    continue;
                }
                else if (v_neighbors[i_idx] == onetight_list[j_idx]) {
                    i_idx++;
                    j_idx++;
                    continue;
                }

                // if this point is reached, this means we found the pair {v, w}
                // we were looking for !!
                int w = onetight_list[j_idx];

                int weight_v = (*weights)[v];
                int weight_w = (*weights)[w];

                if (x_weight >= weight_v + weight_w) {
                    i_idx++;
                    continue;
                }

                removeVertex(x);
                addVertex(v);
                addVertex(w);
                return true;
            }
        } // for(int v : onetight_list) {
    } // for(int x : cadidate_list) {

    return false;
} // bool SolutionMIS::twoImprovment()

bool SolutionMIS::threeImprovement()
{
    assert(isMaximal());

    // for each 2-tight vertex u..
    for (int idx = solution_size_; idx < g_->GetSize(); idx++) {
        int u = solution_[idx];
        if (tightness_[u] != 2) continue;

        // temporarly remove neighbors vertices x and y from the solution
        vector<int> xy;
        for (int j : g_->GetAdj(u)) {
            if (position_[j] < solution_size_) {
                xy.push_back(j);
            }
        }
        removeVertex(xy[0]);
        removeVertex(xy[1]);

        int weight_xy = (*weights)[xy[0]] + (*weights)[xy[1]];

        // temporarly add vertex u (which now is free)
        addVertex(u);

        // if there are less than two free vertices, we are done with u
        if (free_size_ >= 2) {
            // temporarly add each free vertex that is neighbor of x
            for (int v : g_->GetAdj(xy[0])) {
                if (position_[v] >= solution_size_ && position_[v] < solution_size_ + free_size_) {
                    addVertex(v);
                    // if the solution is not maximal, adding any free vertex w will
                    // create a valid solution (thus leading to a 3-improvement)
                    if (!isMaximal()) {
                        int weight_uvz = (*weights)[u] + (*weights)[v] +
                            (*weights)[solution_[solution_size_]];
                        if (weight_uvz > weight_xy) {
                            addVertex(solution_[solution_size_]);
                            return true;
                        }
                    }
                    // remove back v
                    removeVertex(v);
                }
            }
        }

        // add back x, and y, and remove u
        removeVertex(u);
        addVertex(xy[0]);
        addVertex(xy[1]);
    } // for (size_t idx = solution_size_; idx < g_->GetSize(); idx++)

    return false;
} // bool SolutionMIS::threeImprovement()

void SolutionMIS::force(int k)
{
    for (int i = 0; i < k; i++) {
        // select a non solution vertex to add
        int nonsolution_size = g_->GetSize() - (solution_size_ + free_size_);
        uniform_int_distribution<int> discrete_distribution(0, nonsolution_size - 1);
        int nonsolution_pos = discrete_distribution(generator);
        int vertex = solution_[solution_size_ + free_size_ + nonsolution_pos];

        // remove the neighboring vertices as necessary
        for (int neighbor : g_->GetAdj(vertex)) {
            if (position_[neighbor] < solution_size_) {
                removeVertex(neighbor);
            }
        }
        addVertex(vertex);
    }
} // void SolutionMIS::force()
