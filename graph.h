#ifndef GRAPH_H
#define GRAPH_H

#include <queue>
#include <list>
#include <cstdlib>
#include <cmath>

#include "utils/rational.h"

using namespace std;    

class static_flow_graph {
private:	
	int *head, *next, *to, *prev, *h, *cur_edge, *active, *inactive, *prevB, *nextB;
	double *cap, *cap_a, *cap_b, *flow, *ex;
	bool *g;
	bool reversed;
	int edge, n, m, s, t;
	int maxheight;
	double eps;
	double lambda;
	queue<int> active_vrtx;
	long long pushCounter;
	long long relabelCounter;
	long long globalUpdates;
	long long gapCounter;
	long long workDone;
	long long globalUpdateBarrier;
	long long nextUpdate;

	/**
	 * Static constants.
	 */
	static const int GLOBAL_RELABEL_ARC_COEFF = 1;
	static const int GLOBAL_RELABEL_NODE_COEFF = 6;

	static const int RELABEL_WORK_CONST = 12;
	static const int RELABEL_WORK_PER_ARC = 1;

	static const int GLOBAL_UPDATE_WORK_CONST   = 5;
	static const int GLOBAL_UPDATE_WORK_PER_ARC = 2;

	static const int PUSH_WORK_CONST = 0;

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
	 * Updates the global parameter of capacity functions.
	 *
	 * Essentially the same as above one. This one should be prefferable.
	 */
	void update_lambda(double lambda);
	
	/**
	 * Pushes a flow on arc 'e' incident to 'v' (see preflow-push algorithm).
	 */
	void push(int v, int e);

	/**
	 * Assigns a new height label to a vertex (see preflow-push algorithm).
	 */
	void relabel(int v, int stage);

    /**
     * Gets the node for the next discharge.
     */
    int node_for_discharge(int stage);

	/**
	 * Discharges a vertex (see preflow-push algorithm)
	 */
	void discharge(int v, int stage);	
	
	/**
	 * Adds a single arc (a, b) without a reverse arc.
	 */
	void _add_edge(int a, int b, double A, double B);

	/**
	 * Performs gap heuristic.
	 */	
	void gap(int height);
	/**
	 * Performs the global relabeling heuristic.
	 */
	 void global_relabeling(int stage);

public:
	/**
	 * Initializes graph with
	 * N vertices,
	 * M edges / arcs,
	 * source node S,
	 * sink vertex T
	 * and feasible error e.
	 */
	void initialize_graph(int N, int M, int S, int T);//, double e);

	/**
	 * Adds an arc (a, b) with linear depending capacity function of form At+B.
	 */	
	void add_edge(int a, int b, double A, double B);
	
	/**
	 * Discharges the vertices.
	 */    
    double max_flow();

    double max_flow(int stage);

	
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

	void print_heights(int priority);

	/**
	 * Outputs stats of last flow calculation.
	 */
	void print_stats(int priority);

	/**
	 * Returns the number of pushes last comptutation of flow has used.
	 */
	long long pushes();

	/**
	 * Returns the number of relabels last computation of flow has used.
	 */
	long long relabels();

	/**
	 * Returns the number of gap heuristics used in last compuation of flow.
	 */
	long long gaps();

	/**
	 * Returns the nubmer of global relabeling heuristic used for last
	 * computation of flow.
	 */
	long long global_updates();
	
	void rem(int &head, int node);

	void add_front(int &head, int node);
};

#endif