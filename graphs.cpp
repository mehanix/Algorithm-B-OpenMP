#include "graphs.hpp"
#include <iostream>
#include <queue>
//constructors

Graph::Graph(bool d,bool w) : is_directed(d), is_weighted(w), m(0) {}

Graph::Graph() : is_directed(0), is_weighted(0), m(0) {}

Graph::Graph(const Graph& G) : is_directed(G.directed()), is_weighted(G.weighted()), m(0) {
	for(Vertex v : G.vertex_list()) add_node(v);
	for(auto e : G.edges()) add_edge(e.first,e.second,G.weight(e.first,e.second));
}

//edge access

list< pair<Vertex, Vertex> > Graph::edges(void) const {
	list< pair<Vertex,Vertex> > _all;
	for(auto Adj: Out){
		Vertex v = Adj.first; //non-isolated vertices
		for(Vertex u: Adj.second._list){ //neighbours, in order
			if(directed() || v < u)
			 _all.push_back(pair<Vertex,Vertex>(v,u));
		}
	}
	return _all;
}

bool Graph::adjacent(Vertex u, Vertex v) const {
	if(Out.find(u) == Out.end())
	 return 0;
	return (Out.at(u)._map.find(v) != Out.at(u)._map.end());
}

double Graph::weight(Vertex u,Vertex v) const {
	if(!adjacent(u,v)) return numeric_limits<double>::infinity();
	return Out.at(u)._map.at(v)._weight;
}

//vertex access

list<Vertex> Graph::vertex_list(void) const {
	list<Vertex> _all;
	for(Vertex v: vertices)
	 _all.push_back(v);
	return _all;
}

int Graph::degree(Vertex v) const {
	if(vertices.find(v) == vertices.end())
	 return -1;
	else if(directed()) return out_degree(v)+in_degree(v);
	else return out_degree(v);
}

int Graph::out_degree(Vertex v) const{
	if(vertices.find(v) == vertices.end()) return -1;
	else if(Out.find(v) == Out.end()) return 0;
	else return Out.at(v)._list.size();
}

int Graph::in_degree(Vertex v) const{
	if(!directed()) return out_degree(v);
	else if(vertices.find(v) == vertices.end()) return -1;
	else if(In.find(v) == In.end()) return 0;
	else return In.at(v)._list.size();
}

list<Vertex> Graph::neighbours(Vertex v) const {
	list<Vertex> out_list;
	if(Out.find(v) != Out.end())
	 for(Vertex u : Out.at(v)._list)
	  out_list.push_back(u);
	return out_list;
}

list<Vertex> Graph::in_neighbours(Vertex v) const {
	if(!directed()) return neighbours(v);
	list<Vertex> in_list;
	if(In.find(v) != In.end())
	 for(Vertex u : In.at(v)._list)
	  in_list.push_back(u);
	return in_list;
}

//dynamic operations

void Graph::insert(Graph::AdjList& Nu, Vertex v, double w){
	Nu._list.push_front(v);
	Graph::Edge e; e._weight = w; e._pos = Nu._list.begin();
	Nu._map[v] = e;
}

bool Graph::add_edge(Vertex u,Vertex v, double w){
	if(vertices.find(u) == vertices.end() || vertices.find(v) == vertices.end() || adjacent(u,v) || u==v)
	 return 0;
	if(!weighted() && w != 1) {
        cout<<u<<' '<<v<<' '<<w<<'\n';
        cout<<"error:not weighted but weight != 1\n";
        return 0;
	}

	insert(Out[u],v,w);

	if(directed())
	 insert(In[v],u,w);
	else
	 insert(Out[v],u,w);

	m++;
	return 1;
}

void Graph::erase(Graph::AdjList& Nu, Vertex v){
	Graph::Edge e = Nu._map[v];
	Nu._map.erase(v);
	Nu._list.erase(e._pos);
}

bool Graph::remove_edge(Vertex u, Vertex v) {
	if(!adjacent(u,v)) return 0;

	erase(Out[u],v);
	if(out_degree(u)==0) Out.erase(u);

	if(directed()){
	 erase(In[v],u);
	 if(in_degree(v)==0) In.erase(v);
	} else{
	 erase(Out[v],u);
	 if(out_degree(v)==0) Out.erase(v);
	}

	m--;
	return 1;
}

bool Graph::add_node(Vertex v) {
	if(vertices.find(v) != vertices.end()) return 0;
	vertices.insert(v);
	return 1;
}

bool Graph::remove_node(Vertex v) {
	if(vertices.find(v) == vertices.end()) return 0;
	for(Vertex u: neighbours(v))
	 remove_edge(v,u);
	for(Vertex u: in_neighbours(v))
	 remove_edge(u,v);
	vertices.erase(v);
	dead_vertices.push_back(v);
	return 1;
}

void Graph::make_first(Graph::AdjList& Nu, Vertex v){
	Nu._list.push_front(v);
	Nu._list.erase(Nu._map[v]._pos);
	Nu._map[v]._pos = Nu._list.begin();
}

void Graph::reorder(list<Vertex> perm){
	perm.reverse();
	for(Vertex v : perm){

		for(Vertex u: in_neighbours(v))
		 make_first(Out[u],v);

		if(directed())
		 for(Vertex u: neighbours(v))
		  make_first(In[u],v);

	}
}

bool Graph::contract(Vertex x, Vertex y, Vertex z) {

	if(!adjacent(x,y)) return 0;
	if(z != x && z != y && vertices.find(z) != vertices.end()) return 0;

	list<Vertex> out_x = neighbours(x), in_x = in_neighbours(x),
					   out_y = neighbours(y), in_y = in_neighbours(y);
	remove_node(x); remove_node(y);
	add_node(z);
	for(Vertex v: out_x) add_edge(z,v);
	for(Vertex v: out_y) add_edge(z,v);
	for(Vertex u: in_x) add_edge(u,z);
	for(Vertex u: in_y) add_edge(u,z);
	return z;
}

int Graph::contract(Vertex u, Vertex v){
	Vertex z = (degree(u) <= degree(v)) ? u : v;
	if(contract(u,v,z)) return z;
	return -1;
}

Vertex Graph::gen_id(void) {
	if(dead_vertices.empty()) return order();
	Vertex zombie = dead_vertices.front();
	dead_vertices.pop_front();
	return zombie;
}

Vertex Graph::add_node(void){
	Vertex v = gen_id();
	add_node(v);
	return v;
}

unsigned int Graph::heuristic(const Graph& graph, Vertex node, Vertex goal) {
    return 0;
}

void Graph::AStarSearch(const Graph& graph, Vertex start, Vertex goal) {
    // Define a struct for storing the information of a node in the search
    struct Node {
        Vertex vertex;
        double g_score;
        double f_score;

        Node(Vertex v, double g, double f) : vertex(v), g_score(g), f_score(f) {}

        // Define operator< for the priority queue
        bool operator<(const Node& other) const {
            return f_score > other.f_score;  // Min-heap based on f_score
        }
    };

    // Priority queue for storing the nodes to be explored
    std::priority_queue<Node> open_set;

    // Hash table for tracking the parent nodes
    std::unordered_map<Vertex, Vertex> parents;

    // Hash table for tracking the g_scores
    std::unordered_map<Vertex, double> g_scores;

    // Hash table for tracking the f_scores
    std::unordered_map<Vertex, double> f_scores;

    // Initialize the starting node
    g_scores[start] = 0.0;
    f_scores[start] = Graph::heuristic(graph, start, goal);
    open_set.emplace(start, 0.0, f_scores[start]);

    while (!open_set.empty()) {
        // Get the node with the minimum f_score
        Node current = open_set.top();
        open_set.pop();

        if (current.vertex == goal) {
            // Goal reached, reconstruct and print the path
            std::vector<Vertex> path;
            Vertex v = goal;
            while (v != start) {
                path.push_back(v);
                v = parents[v];
            }
            path.push_back(start);

            std::cout << "Path: ";
            for (auto it = path.rbegin(); it != path.rend(); ++it)
                std::cout << *it << " ";
            std::cout << std::endl;

            return;
        }

        // Explore the neighbors of the current node
        for (Vertex neighbor : graph.neighbours(current.vertex)) {
            double tentative_g_score = g_scores[current.vertex] + graph.weight(current.vertex, neighbor);

            // Check if the neighbor has not been evaluated or the new path is better
            if (g_scores.find(neighbor) == g_scores.end() || tentative_g_score < g_scores[neighbor]) {
                parents[neighbor] = current.vertex;
                g_scores[neighbor] = tentative_g_score;
                f_scores[neighbor] = tentative_g_score + Graph::heuristic(graph, neighbor, goal);
                open_set.emplace(neighbor, tentative_g_score, f_scores[neighbor]);
            }
        }
    }

    // No path found
    std::cout << "No path found from " << start << " to " << goal << std::endl;
}

#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>
#include <limits>

// Function to perform Dijkstra's algorithm on a graph
void Graph::Dijkstra(const Graph& graph, Vertex start) {
    // Define a struct for storing the information of a node in the search
    struct Node {
        Vertex vertex;
        double distance;

        Node(Vertex v, double d) : vertex(v), distance(d) {}

        // Define operator< for the priority queue
        bool operator<(const Node& other) const {
            return distance > other.distance;  // Min-heap based on distance
        }
    };

    // Priority queue for storing the nodes to be explored
    std::priority_queue<Node> open_set;

    // Hash table for tracking the parent nodes
    std::unordered_map<Vertex, Vertex> parents;

    // Hash table for tracking the distances
    std::unordered_map<Vertex, double> distances;

    // Initialize the distances
    for (Vertex v : graph.vertex_list()) {
        if (v == start)
            distances[v] = 0.0;
        else
            distances[v] = std::numeric_limits<double>::infinity();
    }

    open_set.emplace(start, 0.0);

    while (!open_set.empty()) {
        // Get the node with the minimum distance
        Node current = open_set.top();
        open_set.pop();

        // Explore the neighbors of the current node
        for (Vertex neighbor : graph.neighbours(current.vertex)) {
            double distance = current.distance + graph.weight(current.vertex, neighbor);

            // Update the distance if a shorter path is found
            if (distance < distances[neighbor]) {
                parents[neighbor] = current.vertex;
                distances[neighbor] = distance;
                open_set.emplace(neighbor, distance);
            }
        }
    }

    // Print the shortest distances and paths from the start vertex
    for (Vertex v : graph.vertex_list()) {
        std::cout << "Shortest distance from " << start << " to " << v << ": ";
        if (distances[v] == std::numeric_limits<double>::infinity()) {
            std::cout << "No path found" << std::endl;
        } else {
            std::cout << distances[v] << std::endl;

            // Reconstruct and print the path
            std::vector<Vertex> path;
            Vertex node = v;
            while (node != start) {
                path.push_back(node);
                node = parents[node];
            }
            path.push_back(start);

            std::cout << "Path: ";
            for (auto it = path.rbegin(); it != path.rend(); ++it)
                std::cout << *it << " ";
            std::cout << std::endl;
        }
    }
}

