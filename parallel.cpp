#include "graphs.hpp"
#include "bag.hpp"
#include <omp.h>
#include <iostream>
#include <algorithm> 
#include <fstream>
#include <cstring>
using namespace std;

/** Parallel implementation of turning undirected graph G to a directed graph, from class. */
void ParallelGraphAlgorithms::make_directed(Graph &G)
{
	if (!G.directed())
	{
		// set flag
		G.is_directed = 1;

		// as the graph is now directed, we now use the In-neighbourhood adjacency list, so we create it.
		for (unsigned int v : G.vertex_list())
			G.In[v]; // create inlist

		// for each vertex v
		for (unsigned int v : G.vertex_list())
		{
			if (G.degree(v) > 0)
			{
				// for each neighbour of v, in parallel, get a reference to the neighbour
				// and add the opposing edge in the in-neighbourhood
				#pragma omp parallel for
				for (unsigned int i = 0; i < G.Out[v]._map.bucket_count(); i++)
				{
					for (auto it = G.Out[v]._map.begin(i); it != G.Out[v]._map.end(i); it++)
					{
						#pragma omp critical
						G.insert(G.In[it->first], v, 1);
					}
				}
			}
		}
	}
}

void ParallelGraphAlgorithms::algorithm_B(Graph G, vector<unsigned int> &p)
{
	/** Step 1. Use the make_directed function from class to turn G into a weighted graph.*/
	make_directed(G);

	/** Step 1.1. Initialise parent array according to the LU framework definition.
	 * Each vertex is initially its own parent. */
	#pragma omp parallel for
	for (unsigned int i = 0; i < G.Out.bucket_count(); i++)
	{
		for (auto it = G.Out.begin(i); it != G.Out.end(i); it++)
		{
			unsigned int v = it->first;
			p[v] = v;
		}
	}

	/** We loop continuously until the end condition has been met.*/
	bool done(false);
	while (!done)
	{
		done = true;

		/** Step 2.1. Consider all arcs w->w in parallel **/
		#pragma omp parallel for
		for (unsigned int i = 0; i < G.Out.bucket_count(); i++)
		{
			for (auto itv = G.Out.begin(i); itv != G.Out.end(i); itv++)
			{
				/** get v's value*/
				unsigned int v = itv->first;

				/** determine v's parent using min rule with reducer*/
				unsigned int parent = p[v];
				#pragma omp parallel for reduction(min : parent)
				for (unsigned int j = 0; j < G.Out[v]._map.bucket_count(); j++)
				{
					for (auto ite = G.Out[v]._map.begin(j); ite != G.Out[v]._map.end(j); ite++)
					{
						/** get w's value */
						unsigned int w = ite->first;
						if (v > w)
						{
							/** Set the parent using the min rule*/
							parent = min(parent, w);
						}
					}
				}
				/** Finally, v's parent in p[v]*/
				p[v] = parent;
			}
		}

		/** Step 3.
		 * Since adding/removing edges from a graph while iterating through said edges is a bad idea
		 * (leads to losing pointer reference/ends up iterating through more edges than it should)
		 * I do this step in two passes:
		 * 	3.0: Create 2 bitmaps where I will flag which edges should be added/deleted.
		 *  3.1: For each edge in parallel, flag for deletion if w=v.p. Else, flag for adding to the graph.
		 *  3.2: Iterate through the bitmaps in parallel and do the deletion/adding of edges without worrying about losing
		 *       pointer references/iterating through recently added edges. 
		 *  
		*/

		/** Step 3.0: */
		unsigned int n = G.order();
		bool bitmap_todelete[n][n];
		bool bitmap_toadd[n][n];
		memset(bitmap_todelete, 0, n * n * sizeof(bool));
		memset(bitmap_toadd, 0, n * n * sizeof(bool));


		/** Step 3.1: */
		// for each arc v->w in parallel
		#pragma omp parallel for
		for (unsigned int i = 0; i < G.Out.bucket_count(); i++)
		{
			for (auto itv = G.Out.begin(i); itv != G.Out.end(i); itv++)
			{
				unsigned int v = itv->first;
				#pragma omp parallel for
				for (unsigned int j = 0; j < G.Out[v]._map.bucket_count(); j++)
				{
					for (auto ite = G.Out[v]._map.begin(j); ite != G.Out[v]._map.end(j); ite++)
					{
						unsigned int w = ite->first;
						bitmap_todelete[v][w] = 1; // mark for deletion
						if (p[v] != w)
						{
							// mark for adding. doing this to avoid iterating through the newly added edges in this step.
							bitmap_toadd[w][p[v]] = 1;
						}
					}
				}
			}
		}

		/** Step 3.2: */
		// delete/add marked edges, in parallel
		#pragma omp parallel for
		for (unsigned int i = 0; i < n; i++)
		{
			#pragma omp parallel for
			for (unsigned int j = 0; j < n; j++)
			{
				if (bitmap_todelete[i][j] == 1)
				{
					#pragma omp critical
					G.remove_edge(i, j);
				}
				if (bitmap_toadd[i][j] == 1)
				{
					#pragma omp critical
					G.add_edge(i, j);
				}
			}
		}
		/** Step 4. Finally, we consider all vertices v in parallel. If v!=v.p, add a new arc v->v.p.*/
		#pragma omp parallel for
		for (unsigned int i = 0; i < G.Out.bucket_count(); i++)
		{
			for (auto it = G.Out.begin(i); it != G.Out.end(i); it++)
			{
				unsigned int v = it->first;
				if (p[v] != v)
				{
					#pragma omp critical
					G.add_edge(p[v], v);
				}
			}
		}

		/** Step 5. Check the end condition. For each edge in parallel, if we find one that 
		 * does not have v.p = w.p, with w.p being either v or w, we set the done flag to false,
		 * thus continuing the algorithm.
		*/
		#pragma omp parallel for
		for (unsigned int i = 0; i < G.Out.bucket_count(); i++)
		{
			for (auto itv = G.Out.begin(i); itv != G.Out.end(i); itv++)
			{
				unsigned int v = itv->first;
				#pragma omp parallel for
				for (unsigned int j = 0; j < G.Out[v]._map.bucket_count(); j++)
				{
					for (auto ite = G.Out[v]._map.begin(j); ite != G.Out[v]._map.end(j); ite++)
					{
						unsigned int w = ite->first;
						if (!(p[v] == p[w] && (p[v] == v || p[v] == w)))
						{
							done = false;
						}
					}
				}
			}
		}
	}

	/** We're done! */
}
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// Algorithm S done in the lab.
// Used for checking algorithm B's output against a source of truth.
bool ParallelGraphAlgorithms::shortcut(Graph& G, vector<unsigned int>& p, vector<unsigned int>& o){

	bool changes(0);

	#pragma omp parallel for
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){
		for(auto it = G.vertices.begin(i); it != G.vertices.end(i); it++){
			unsigned int v = *it;
			o[v] = p[v];
		}
	}

	#pragma omp parallel for
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){
		for(auto it = G.vertices.begin(i); it != G.vertices.end(i); it++){
			unsigned int v = *it;
			p[v] = o[o[v]];
			if(p[v] != o[v]) changes = 1;
		}
	}

	return changes;
}


bool ParallelGraphAlgorithms::parent_connect(Graph& G, vector<unsigned int>& p, vector<unsigned int>& o){

	bool changes(0);

	/* CMT: atribuire in paralel */

	#pragma omp parallel for
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){
		for(auto itv = G.vertices.begin(i); itv != G.vertices.end(i); itv++){
			unsigned int v = *itv;
			o[v] = p[v];
		}
	}

	#pragma omp parallel for
	for(unsigned int i = 0; i < G.Out.bucket_count(); i++){
		for(auto itv = G.Out.begin(i); itv != G.Out.end(i); itv++){

			unsigned int v = itv->first;

			#pragma omp parallel for
			for(unsigned int j = 0; j < G.Out[v]._map.bucket_count(); j++){
				for(auto ite = G.Out[v]._map.begin(j); ite != G.Out[v]._map.end(j); ite++){

					unsigned int w = ite->first;

					if(o[v] > o[w])
					 #pragma omp critical
					 p[o[v]] = min(p[o[v]],o[w]);

				}
			}

		}
	}

	return changes;
}

void ParallelGraphAlgorithms::algorithm_S(Graph G, vector<unsigned int>& p, vector<unsigned int>& o){

	#pragma omp parallel for
	for(unsigned int i = 0; i < G.Out.bucket_count(); i++){
		for(auto it = G.Out.begin(i); it != G.Out.end(i); it++){
			unsigned int v = it->first;
			p[v] = v;
		}
	}

	bool not_done(1);
	while(not_done){
		not_done = parent_connect(G,p,o);
		bool not_star(1);
		while(not_star) not_star = shortcut(G,p,o);
	}

}