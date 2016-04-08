#ifndef GRAPH_H
#define GRAPH_H

#include <queue>
#include <cstdlib>
#include <cmath>

using namespace std;    

class static_flow_graph {
private:	
	int *head, *next, *to, *prev, *h, *cur_edge;
	double *cap, *cap_a, *cap_b, *flow, *ex;
	bool *g;
	bool reversed;
	int edge, n, m, s, t;
	double eps;
	double lambda;
	queue<int> active_vrtx;
	int pushCounter;
	int relabelCounter;

	/**
	 * Initializes preflow-push algorithm
	 */
	void initialize_flow(double L);
	

	/**
	 * Decreases the global parameter of capacity functions.
	 *
	 * lambda_dif should be non-positive.
	 */
	void decrease_cap(double lambda_dif);
	
	/**
	 * Pushes a flow on arc 'e' incident to 'v' (see preflow-push algorithm).
	 */
	void push(int v, int e, int s, int t);

	/**
	 * Assigns a new height label to a vertex (see preflow-push algorithm).
	 */
	void relabel(int v);                              	
	
	/**
	 * Discharges a vertex (see preflow-push algorithm)
	 */
	void discharge(int v, int s, int t);	
	
	/**
	 * Adds a single arc (a, b) without a reverse arc.
	 */
	void _add_edge(int a, int b, double A, double B);


public:
	/**
	 * Initializes graph with
	 * N vertices,
	 * M edges / arcs,
	 * source node S,
	 * sink vertex T
	 * and feasible error e.
	 */
	void initialize_graph(int N, int M, int S, int T, double e);

	/**
	 * Adds an arc (a, b) with linear depending capacity function of form At+B.
	 */	
	void add_edge(int a, int b, double A, double B);
	
	/**
	 * Discharges the vertices.
	 */    
    double max_flow();

	
	/**
	 * Reverses current flow state. That is, edges and their capacity and flow are reoriented 
	 * backwards, sink and source nodes are swapped.
	 */
	void reverse_flow();

	/**
	 * Returns true if one can send additional flor from parameter vertex v to sink vertex t.
	 */	
	bool has_flow(int v);
	
	/**
	 * Finds the leftmost breakpoint for a monotomic flow (see GGT method for parametric 
	 * flows for details).
	 * Note that leftmost breakpoint for monotomic flow corresponds to cut consisting of 
	 * single node {t}.
	 *
	 * Also, init_lambad is required to be positive non-zero and greater then leftmost 
	 * breakpoint. Warning!!! Function might not work as expected with init_lambda 
	 * values close to zero (the closer the value to zero the more impact calculation 
	 * errors will have).
	 */
	double leftmost_breakpoint(double init_lambda);

	double get_flow(int idx);

	void show(int priority);

	/**
	 * Returns the number of pushes last comptutation of flow has used.
	 */
	int pushes();

	/**
	 * Returns the number of relabels last computation of flow has used.
	 */
	int relabels();
};

#endif