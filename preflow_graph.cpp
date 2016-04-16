#include "graph.h"
#include "utils/logger.h"
#include <math.h>

void static_flow_graph::initialize_flow(double L) {
	pushCounter = 0;
	relabelCounter = 0;
	workDone = 0;
	globalUpdates = 0;
	lambda = L;
	double x;
	for(int i = 0; i < m; ++i) {
		cap[i] = cap_a[i] + lambda * cap_b[i];
		flow[i] = 0;
	}
	globalUpdateBarrier = (int) (GLOBAL_RELABEL_ARC_COEFF * m + GLOBAL_RELABEL_NODE_COEFF * n);
}

/**
 * Breaks active vertices structure.
 * Have to be combined with global_relabeling.
 */
void static_flow_graph::decrease_cap(double lambda_dif) {
	for(int v = 0; v < n; ++v) {
		for(int i = head[v]; i != -1; i = next[i]) {
			cap[i] += lambda_dif * cap_b[i];
			if(flow[i] <= cap[i]) continue;
			double val = flow[i] - cap[i];

            ex[v] += val;
			ex[to[i]] -= val;
			flow[i] -= val;
			flow[i^1] = -flow[i];
		}
	}
}

/**
 * Breaks active vertices structure.
 * Have to be combined with global_relabeling.
 */
void static_flow_graph::update_lambda(double lambda) {
	for(int v = 0; v < n; ++v) {
		ex[v] = 0;
		for(int i = head[v]; i != -1; i = next[i]) {
			cap[i] = cap_a[i] + lambda * cap_b[i];
			if(flow[i] > cap[i]) {
				flow[i] = cap[i];
				flow[i^1] = -cap[i];
			}
			ex[v] -= flow[i];			
		}
	}
}

void static_flow_graph::push(int v, int e, int s, int t) {
	pushCounter++;
	workDone += PUSH_WORK_CONST;
#ifdef DEBUG
	print("Pushing from ", 3);
	print(v, 3);
	print(" to ", 3);
	print(to[e], 3);
	print("\n", 3);
#endif
	if(to[e] != s && to[e] != t && ex[to[e]] < eps) {
		rem(inactive[h[to[e]]], to[e]);
		maxheight = max(maxheight, h[to[e]]);
		add_front(active[h[to[e]]], to[e]);
	}
	double val = (ex[v] < cap[e] - flow[e] ? ex[v] : cap[e] - flow[e] );
    ex[v] -= val;
	ex[to[e]] += val;
	flow[e] += val;
	flow[e^1] = -flow[e];
	if(ex[v] < eps) {
		rem(active[h[v]], v);
		add_front(inactive[h[v]], v);
	}	
}

void static_flow_graph::global_relabeling() {
#ifdef DEBUG	
	print("Starting global relabeling\n", 1);
#endif
	globalUpdates++;
	queue<int> q;
	q.push(t);
	for(int i = 0; i < n; ++i) {
		g[i] = (i == t) || (i == s);
	}
	for(int i = 0; i < 2 * n; ++i) {
		active[i] = inactive[i] = -1;
	}
	maxheight = 0;
	while(!q.empty()) {
		int temp = q.front();
		q.pop();
		for(int i = head[temp]; i != -1; i = next[i]) {
			if(!g[to[i]] && (cap[i^1] - flow[i^1] > eps)) {
				h[to[i]] = h[temp] + 1;
				g[to[i]] = 1;
				q.push(to[i]);
				if(ex[to[i]] > 0) {
					add_front(active[h[to[i]]], to[i]);
					maxheight = h[to[i]];
				} else {
					add_front(inactive[h[to[i]]], to[i]);
				}
			}
		}
	}
	q.push(s);
	while(!q.empty()) {
		int temp = q.front();
		q.pop();
		for(int i = head[temp]; i != -1; i = next[i]) {
			if(!g[to[i]] && (cap[i^1] - flow[i^1] > eps)) {
				h[to[i]] = h[temp] + 1;
				g[to[i]] = 1;
				if(ex[to[i]] > 0) {
					add_front(active[h[to[i]]], to[i]);
					maxheight = h[to[i]];
				} else {
					add_front(inactive[h[to[i]]], to[i]);
				}
				q.push(to[i]);
			}
		}
	}
}

void static_flow_graph::relabel(int v) {                              	
#ifdef DEBUG	
	print("Relabeling ", 3);
	print(v, 3);
#endif	
	relabelCounter++;
	workDone += RELABEL_WORK_CONST;
	int u = -1;
	for(int i = head[v]; i != -1; i = next[i]) {	
		if(cap[i] - flow[i] > eps && (u == -1 || h[to[u]] > h[to[i]])) {
			u = i;
		}
		workDone += RELABEL_WORK_PER_ARC;
	}
#ifdef DEBUG
	print("\nNew label is ", 3);
	print(h[to[u]] + 1, 3);
	print("\n", 3);
#endif	
	rem(active[h[v]], v);
	cur_edge[v] = u;
	h[v] = h[to[u]] + 1;	
	add_front(active[h[v]], v);
	maxheight = max(maxheight, h[v]);
}
	
void static_flow_graph::discharge(int v, int s, int t) {
	if(head[v] == -1) return;
	while(1) {
		if(cap[cur_edge[v]] - flow[cur_edge[v]] > eps && h[v] == h[to[cur_edge[v]]] + 1) {
			push(v, cur_edge[v], s, t);
		}
		if(ex[v] == 0) break;
		cur_edge[v] = next[cur_edge[v]];
		if(cur_edge[v] == -1) {
			relabel(v);
			cur_edge[v] = head[v];
		}
	}
}
	
void static_flow_graph::_add_edge(int a, int b, double A, double B) {
	if(edge == m) {
		//Insert error message
		return;
	}
	next[edge] = head[a];
	head[a] = edge;
	to[edge] = b;
	cap_a[edge] = A;
	cap_b[edge] = B;
	edge++;
}


void static_flow_graph::initialize_graph(int N, int M, int S, int T) {
   	head = new int[N];
   	prev = new int[N];
   	nextB = new int[N];
   	prevB = new int[N];
   	active = new int [2 * N];
   	inactive = new int [2 * N];
   	g = new bool[N];
   	cur_edge = new int[N];
   	ex = new double[N];
   	h = new int[N];

    next = new int[2 * M];
   	to = new int[2 * M];
    for(int i = 0; i < N; ++i) {
    	head[i] = -1;
    }
    cap_a = new double[2 * M];
    cap_b = new double[2 * M];
    cap = new double[2 * M];
    flow = new double[2 * M];
        
    reversed = false;
    n = N;
    m = 2 * M;
    s = S;
    t = T;
    edge = 0;    			
	//OBSOLETTE
    eps = 1e-9;
}

void static_flow_graph::add_edge(int a, int b, double A, double B) {
	_add_edge(a, b, A, B);
	_add_edge(b, a, 0, 0);		
}     
    
double static_flow_graph::max_flow() {
#ifdef DEBUG
	print("Starting maxflow\n", 3);
#endif	
	while(1) {
		while(maxheight > -1 && active[maxheight] == -1) {
			maxheight--;
		}
		if(maxheight == -1) break;
#ifdef DEBUG		
		print("Max height = ", 3);
		print(maxheight, 3);
		print("\n", 3);
#endif
		discharge(active[maxheight], s, t);
		if(workDone >= nextUpdate) {
			global_relabeling();
			nextUpdate += globalUpdateBarrier;
		}
		//print_heights(4);
	}
	return ex[t];	
}

	
/**
 * Inverses current flow state. That is, edges and their capacity and flow are reoriented 
 * backwards, sink and source nodes are swapped.
 */
void static_flow_graph::reverse_flow() {
	reversed = !reversed;
	swap(s, t);
	for(int i = 0; i < m / 2; ++i) {
		swap(cap[2 * i], cap[2 * i + 1]);
		swap(cap_a[2 * i], cap_a[2 * i + 1]);
        swap(cap_b[2 * i], cap_b[2 * i + 1]);
	    swap(flow[2 * i], flow[2 * i + 1]);
	} 
}
	
bool static_flow_graph::has_flow(int v) {
	if(v == t) return true;
	if(g[v]) return false;
	g[v] = true;

	for(int i = head[v]; i != -1; i = next[i]) {
		if(cap[i] - flow[i] > eps && has_flow(to[i])) return true;
	}
	return false;
}
	
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
double static_flow_graph::leftmost_breakpoint(double init_lambda) {
	initialize_flow(init_lambda);
	reverse_flow();
	
	for(int i = 0; i < n; ++i) {
		h[i] = n * (i == s);
		ex[i] = 0;
		cur_edge[i] = head[i];
	}
	for(int i = head[s]; i != -1; i = next[i]) {
		ex[to[i]] = cap[i];
		flow[i] = cap[i];
		flow[i^1] = -flow[i];
	}
	maxheight = 0;
	lambda = init_lambda;
	global_relabeling();
	nextUpdate += globalUpdateBarrier;
		
    double fl = max_flow();
    double load = 0;
    for(int i = head[s]; i != -1; i = next[i]) {
    	load += cap[i];
    }

	while(fl - load < eps) {
		queue<int> q;
		reverse_flow();
		for(int i = 0; i < n; i++) {
			g[i] = (i == s);
		}			
		q.push(s);
		
		int cut_size = 0;
		
		while(!q.empty()) {
			cut_size++;
			int temp = q.front();
			//cerr << "q:" << temp << endl;
			q.pop();
			for(int i = head[temp]; i != -1; i = next[i]) {
				if(cap[i] - flow[i] > eps && !g[to[i]]) {
					q.push(to[i]);
					g[to[i]] = true;
				}
			}
		}
		/**
		 * PARANOIC
		 */
		if(cut_size == 1) break;
	    g[s] = false;
		double mult = 0, stat = 0;
		for(int i = 0; i < n; ++i) {
			if(g[i]) {
				for(int j = head[i]; j != -1; j = next[j]) {
					if(!g[to[j]]) {
						stat -= cap_a[j];
						mult += cap_b[j];
					}
				}
			}
			if(i == s) {
				for(int j = head[i]; j != -1; j = next[j]) {
					if(g[to[j]]) {
						stat += cap_a[j];
						mult -= cap_b[j];
					}
				}
			}
		}
		reverse_flow();

#ifdef DEBUG
		print(stat, 1);
		print(" ", 1);
		print(mult, 1);
		print(" New L = ", 1);
#endif
		double new_lambda = stat / mult;
		//cerr << s_load << " " << cut << " " << tau << " " << new_tau << endl;
	    //decrease_cap(new_lambda - lambda);
	    update_lambda(new_lambda);    
	    lambda = new_lambda;
	   	
	   	global_relabeling();
		nextUpdate += globalUpdateBarrier;

	    fl = max_flow();
	    //cerr << 1 / lambda << " " << fl << endl;
	}
	/*
	 * PARANOIC
	 */
	if(reversed) reverse_flow();
	
	/*
	for(int i = 0; i < n; ++i) {
		for(int j = head[i]; j != -1; j = next[j]) {
			cerr << i + 1 << "->" << to[j] + 1 << ": " << flow[j] << " " << cap[j] << endl;
		}
	}
	*/

	return lambda;
}

void static_flow_graph::show(int priority) {
	print("Number of nodes = ", priority);
	print(this->n, priority);
	print("\n", priority);
	print("Number of arcs = ", priority);
	print(m, priority);
	print("\nSource node is ", priority);
	print(s, priority);
	print("\nSink node is ", priority);
	print(t, priority);
	print("\nCurrent lambda = ", priority);
	print(lambda, priority);
	print("\n", priority);
           	
	for(int i = 0; i < n; ++i) {
		print(i, priority);
		print(" ex = ", priority);
		print(ex[i], priority);
		print(":\n", priority);

		for(int j = head[i]; j != -1; j = next[j]) {
			print("arc n ", priority);
			print(j, priority);
			print(" at ", priority);
			print(to[j], priority);
			print(" with cap = ", priority);
			print(cap_a[j], priority);
			print(" + ", priority);
			print(lambda, priority);
			print(" * ", priority);
			print(cap_b[j], priority);
			print(" = ", priority);
			print(cap[j], priority);
			print(" and flow = ", priority);
			print(flow[j], priority);
			print("\n", priority);
			
		}
		print("\n", priority);
	}
}

void static_flow_graph::print_heights(int priority) {
	for(int i = 0; i < 2 * n; ++i) {
		int temp = active[i];
		print("Active vertices at heigth ", priority);
		print(i, priority);
		print(":", priority);
		while(temp != -1) {
			print(temp, priority);
			print(" ", priority);
			temp = nextB[temp];
		}
		print(" (0-indexed)\n", priority);
	}
}

void static_flow_graph::print_stats(int priority) {
	print("Number of pushes: ", priority);
	print(pushCounter, priority);
	print("\n", 2);
	print("Nubmer of relabels: ", priority);
	print(relabelCounter, priority);
	print("\n", 2);
	print("Nubmer of global relabels: ", priority);
	print(globalUpdates, priority);
	print("\n", priority);
}

void static_flow_graph::add_front(int &head, int node) {
	nextB[node] = head;
	prevB[node] = -1;
	if(head != -1) {
		prevB[head] = node;
	}		 
	head = node;
}

void static_flow_graph::rem(int &head, int node) {
	if(prevB[node] != -1) {
		nextB[prevB[node]] = nextB[node];
	}
	if(nextB[node] != -1) {
		prevB[nextB[node]] = prevB[node];
	}
	if(head == node) {
		head = nextB[node];
	}
}                                 

double static_flow_graph::get_flow(int idx) {
	return flow[2 * idx];
}

long long static_flow_graph::pushes() {
	return pushCounter;
}

long long static_flow_graph:: relabels() {
	return relabelCounter;
}

long long static_flow_graph::global_updates() {
	return globalUpdates;
}	
