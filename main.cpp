#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <map>


using namespace std;


struct Point {
    int taxi_id;
    std::string date_time;
    double longitude;
    double latitude;
    int cluster;
    double norm_longitude;
    double norm_latitude;
};

std::vector<Point> read_taxi_data(const std::string &filename) {
    std::ifstream infile(filename);
    std::vector<Point> taxi_data;
    std::string line;
    bool first_line = true;
    while (std::getline(infile, line)) {
        if (first_line) {
            first_line = false;
            continue;
        }
        std::istringstream iss(line);
        Point data;
        std::string field;
        std::getline(iss, field, ',');
        data.taxi_id = std::stoi(field);
        std::getline(iss, data.date_time, ',');
        std::getline(iss, field, ',');
        data.longitude = std::stod(field);
        std::getline(iss, field, ',');
        data.latitude = std::stod(field);
        std::getline(iss, field, ',');
        data.cluster = std::stoi(field);
        std::getline(iss, field, ',');
        data.norm_longitude = std::stod(field);
        std::getline(iss, field, ',');
        data.norm_latitude = std::stod(field);
        taxi_data.push_back(data);
    }
    infile.close();
    return taxi_data;
}


struct GraphEdge {
    int u;
    int v;
    double weight;
};


struct GraphNode {
    int id;
    Eigen::VectorXd pos;
    vector<Point> points;
};

class Graph {
public:
    Graph(const std::map<int, std::vector<Point>> &clusters) {
        // Compute scaler
        std::vector<Point> all_points;
        for (const auto &c: clusters) {
            all_points.insert(all_points.end(), c.second.begin(), c.second.end());
        }
        std::map<std::string, double> scaler = compute_scaler(all_points);

        // Create vertices
        for (const auto &c: clusters) {
            vertices_.push_back(c.first);
        }

        // Compute edge weights
        std::vector<GraphNode> nodes;
        for (const auto &c: clusters) {
            nodes.push_back({c.first, calculate_centroid(c.second), c.second});
        }
        std::vector<GraphEdge> edges = compute_edge_weights(nodes, scaler);

        // Create edges
        for (const auto &e: edges) {
            if (e.weight < 0.3)
                edges_.emplace_back(e.u, e.v);
        }

        vector<int> nbr_edges(vertices().size(), 0);
        for (const auto &item: edges_) {
            nbr_edges[item.first]++;
        }


        vertices_.erase(
                std::remove_if(
                        vertices_.begin(),
                        vertices_.end(),
                        [&nbr_edges, this](const int &n) {
                            if (nbr_edges[n] == 0) {
                                edges_.erase(
                                        std::remove_if(
                                                edges_.begin(),
                                                edges_.end(),
                                                [&n](const std::pair<int, int> &e) {
                                                    return e.first == n || e.second == n;
                                                }),
                                        edges_.end()
                                );
                                return true;
                            }
                            return false;
                        }),
                vertices_.end()
        );

    }

    const std::vector<int> &vertices() const { return vertices_; }

    const std::vector<std::pair<int, int>> &edges() const { return edges_; }

private:
    std::vector<int> vertices_;
    std::vector<std::pair<int, int>> edges_;

    Eigen::VectorXd calculate_centroid(const std::vector<Point> &points) const {
        Eigen::VectorXd centroid(2);
        centroid << 0, 0;
        for (const auto &p: points) {
            centroid[0] += p.norm_longitude;
            centroid[1] += p.norm_latitude;
        }
        centroid /= points.size();
        return centroid;
    }

    std::vector<GraphEdge>
    compute_edge_weights(const std::vector<GraphNode> &nodes, const std::map<std::string, double> &scaler) const {
        std::vector<GraphEdge> edges;
        Eigen::Vector2d u_pos, v_pos;
        double dist;
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                u_pos = nodes[i].pos;
                v_pos = nodes[j].pos;
                dist = (u_pos - v_pos).norm();
                dist = scaler.at("scale") * (dist - scaler.at("mean")) / scaler.at("std");
                edges.push_back({i, j, dist});
            }
        }
        return edges;
    }

    std::map<std::string, double> compute_scaler(const std::vector<Point> &points) const {
        double sum = 0.0, sq_sum = 0.0;
        for (const auto &p: points) {
            sum += p.norm_longitude;
            sq_sum += p.norm_longitude * p.norm_longitude;
            sum += p.norm_latitude;
            sq_sum += p.norm_latitude * p.norm_latitude;
        }
        double mean = sum / (2 * points.size());
        double std = std::sqrt(sq_sum / (2 * points.size()) - mean * mean);
        double scale = 1.0 / std;
        return {{"mean",  mean},
                {"std",   std},
                {"scale", scale}};
    }
};

void save_graph_to_dot(const Graph &graph, const std::string &filename) {
    std::ofstream file(filename);
    file << "graph G {\n";
    for (const auto &e: graph.edges()) {
        file << e.first << " -- " << e.second << ";\n";
    }
    file << "}\n";
    file.close();
}

int main() {
    vector<Point> points = read_taxi_data(
            R"(C:\Users\W0L1D\CLionProjects\ILISI2_atelier_Shortest_path_A_start\saved_clusters_normCoord_45_0.01.csv)");

    map<int, vector<Point>> clusters;
    for (const Point &p: points) {
        if (clusters.find(p.cluster) == clusters.end())
            clusters[p.cluster] = vector<Point>();
        clusters[p.cluster].push_back(p);
    }

    clusters.erase(-1);

    Graph g(clusters);

    cout << "nombre of vertices : " << g.vertices().size() << endl;
    cout << "nombre of edges : " << g.edges().size() << endl;

    vector<int> nbr_edges(g.vertices().size(), 0);
    for (const auto &item: g.edges()) {
        nbr_edges[item.first]++;
    }

    const auto &vertices = (g.vertices());
    const auto &edges = g.edges();

    sort(nbr_edges.begin(), nbr_edges.end());
    double sum = 0;
    for (int i = 0; i < nbr_edges.size(); i++) {
        sum += nbr_edges[i];
    }
    cout << "average nbr of edges : " << sum / nbr_edges.size() << endl;
    for (int i = 0; i < 40; ++i) {
        cout << "min nbr of edges " << i << " : " << nbr_edges[i] << endl;
    }
    cout << "max nbr of edges : " << nbr_edges[nbr_edges.size() - 1] << endl;

    save_graph_to_dot(g, R"(C:\Users\W0L1D\CLionProjects\ILISI2_atelier_Shortest_path_A_start\graph.dot)");
    return 0;
}
