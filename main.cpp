#include <iostream>
#include <fstream>
#include "graphs.hpp"
#include <algorithm>
#include <time.h> /* time */
using namespace std;

int n, m;
Graph g;

/** Reads file input, turns it into graph object */
void read_graph(bool weighted)
{
    string filename;
    int left, right, cost;
    cout << "Enter file name: ";
    cin >> filename;
    ifstream fin(filename);
    fin >> n >> m;

    while (fin >> left >> right)
    {
        g.add_node(left);
        g.add_node(right);

        if (weighted)
        {
            fin >> cost;
            g.add_edge(left, right, cost);
        }
        else
        {
            g.add_edge(left, right);
        }
    }
    fin.close();
    return;
}

int main()
{
    bool is_directed = false;
    bool is_weighted = false;

    /** Reads graph from file according to exam requirements */
    g = Graph(is_directed, is_weighted);
    read_graph(is_weighted);

    /**
     * Step 0.
     * We will use vector p to keep track of each vertex's parent in the connected components algorithm.
     * We need a reference in order to use it in the algorithm later so I declare it here.
     * */
    vector<unsigned int> p(g.order(), 0);

    /**
     * Run Algorithm B on graph g and output the connected components in p.
    */
    ParallelGraphAlgorithms::algorithm_B(g, p);

    /**
     * The connected components are found in the parent array now.
     * Vertex I is part of component rooted at X if p[I]=X.
     * I will create a list of pairs <connected_component_root, vertex>, sort by the connected component's root,
     * and output it to the screen.
     * */
    cout<<"Output of Algorithm B:\n";
    vector<pair<unsigned int,unsigned int>> connected_components;
    for (unsigned int i=0;i<g.order();i++) {
        pair<int, int> a;
        a.first = p[i];
        a.second = i;
        connected_components.push_back(a);
    }
    sort(connected_components.begin(), connected_components.end());
    unsigned int current_component = connected_components[0].first;

    for(auto& component:connected_components) {
        if(component.first == current_component) {
            cout<<component.second<<' ';
        } else {
            cout<<'\n'<<component.second<<' ';
            current_component = component.first;
        }
    }
    cout<<"\n-----------------------------------------------\n";

    cout<<"For reference, here is the output of algoritm S \nfor computing connected components, from the lab:\n";
    vector<unsigned int> pp(g.order(), 0);
    vector<unsigned int> oo(g.order(), 0);

    ParallelGraphAlgorithms::algorithm_S(g, pp, oo);
    vector<pair<unsigned int,unsigned int>> connected_components_S;
    for (unsigned int i=0;i<g.order();i++) {
        pair<unsigned int, unsigned int> a;
        a.first = p[i];
        a.second = i;
        connected_components_S.push_back(a);
    }
    sort(connected_components_S.begin(), connected_components_S.end());
    current_component = connected_components_S[0].first;
    for(auto& component:connected_components) {
        if(component.first == current_component) {
            cout<<component.second<<' ';
        } else {
            cout<<'\n'<<component.second<<' ';
            current_component = component.first;
        }
    }

    cout<<"\n-----------------------------------------------\n";
    return 0;
}
