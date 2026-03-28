#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <vector>
#include <algorithm>
#include <cmath>

// 1D DBSCAN clustering for temporal data
// Finds clusters of temporally close times, returns indices of the densest cluster
// If no cluster found (all noise), returns indices within eps of the earliest time
//
// Args:
//   times: sorted vector of time values
//   eps: maximum temporal distance between adjacent points to be in same cluster (default 50ns)
//   min_pts: minimum number of points to form a valid cluster (default 2)
//
// Returns: vector of indices (into original times array) forming the densest cluster
//          if no cluster found, returns times
//          if times is empty, returns empty vector
inline const std::vector<int>
find_densest_time_cluster(const std::vector<float>& times,
                          float eps = 50.0f,
                          int min_pts = 2)
{
  if (times.empty()) return {};

    int n = times.size();

    // For sorted 1D data, find neighbors by scanning left/right from each point
    auto get_neighbors = [&](int i) -> std::vector<int> {
        std::vector<int> neighbors;
        // Scan left
        int left = i;
        while (left > 0 && times[i] - times[left - 1] <= eps) left--;
        // Scan right
        int right = i;
        while (right < n - 1 && times[right + 1] - times[i] <= eps) right++;
        // Collect range [left, right]
        for (int j = left; j <= right; j++) neighbors.push_back(j);
        return neighbors;
    };

    // -1 = unvisited, -2 = visited but not assigned (noise), >= 0 = cluster id
    std::vector<int> cluster_id(n, -1);
    int cluster_count = 0;
    std::vector<int> cluster_sizes;

    // DBSCAN clustering on sorted 1D data
    for (int i = 0; i < n; i++) {
        if (cluster_id[i] != -1) continue;  // already visited or assigned

        std::vector<int> neighbors = get_neighbors(i);

        if ((int)neighbors.size() < min_pts) {
            // i is not a core point; mark as visited but not assigned
            cluster_id[i] = -2;
            continue;
        }

        // i is a core point: start new cluster
        int current_cluster = cluster_count;
        cluster_count++;
        cluster_sizes.push_back(0);

        // Add core point to cluster
        cluster_id[i] = current_cluster;
        cluster_sizes[current_cluster]++;

        // Expand cluster through core points
        std::vector<int> to_process;
        for (int nb : neighbors) {
          if (cluster_id[nb] == -1 || cluster_id[nb] == -2) {
                to_process.push_back(nb);
            }
        }

        size_t processed_idx = 0;
        while (processed_idx < to_process.size()) {
            int j = to_process[processed_idx];
            processed_idx++;

            if (cluster_id[j] >= 0) continue;  // already assigned to some cluster

            cluster_id[j] = current_cluster;
            cluster_sizes[current_cluster]++;

            // Check if j is also a core point
            std::vector<int> j_neighbors = get_neighbors(j);
            if ((int)j_neighbors.size() >= min_pts) {
                // j is core: add its unvisited neighbors to expansion
                for (int nb : j_neighbors) {
                    if (cluster_id[nb] == -1 || cluster_id[nb] == -2) {
                        to_process.push_back(nb);
                    }
                }
            }
        }
    }

    // Find densest cluster
    if (cluster_count == 0) {
      // No clusters found, return all indices (all points are noise)
      std::vector<int> all_indices(n);
        for (int i = 0; i < n; i++) all_indices[i] = i;
      return all_indices;
    }

    int densest_cluster = 0;
    int max_size = cluster_sizes[0];
    for (int i = 1; i < cluster_count; i++) {
        if (cluster_sizes[i] > max_size) {
            max_size = cluster_sizes[i];
            densest_cluster = i;
        }
    }

    // Collect indices in densest cluster (sorted by time)
    std::vector<int> result;
    for (int i = 0; i < n; i++) {
        if (cluster_id[i] == densest_cluster) {
            result.push_back(i);
        }
    }

    return result;
}

#endif // CLUSTERING_H
