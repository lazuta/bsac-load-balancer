#include "graph.h"
#include "utils/logger.h"
#include <math.h>
//#define DEBUG

const double eps = 1e-12;

bool is_zero(flow_type f) {
#ifdef RATIONAL_FLOW
	return (f == 0);
#endif
#ifdef DOUBLE_FLOW
	return fabs(f) < eps;
#endif
}

void static_flow_graph::initialize_flow(flow_type L) {
	lambda = L;
	for(int i = 0; i < m; ++i) {
		cap[i] = cap_a[i] + lambda * cap_b[i];
		flow[i] = 0;
	}
	for(int i = 0; i < n; ++i) {
		active[i] = inactive[i] = active[i + n] = inactive[i + n] = -1;
		nextB[i] = -1;
		prevB[i] = -1;
		cur_edge[i] = head[i];
	}		
	globalUpdateBarrier = (int) (GLOBAL_RELABEL_ARC_COEFF * m + GLOBAL_RELABEL_NODE_COEFF * n);
}

void static_flow_graph::clear_stats() {
	pushCounter = 0;
	relabelCounter = 0;
	workDone = 0;
	nextUpdate = 0;
	globalUpdates = 0;
	gapCounter = 0;
}

/**
 * Breaks active vertices structure.
 * Have to be combined with global_relabeling.
 */
void static_flow_graph::decrease_cap(flow_type lambda_dif) {
	for(int v = 0; v < n; ++v) {
		for(int i = head[v]; i != -1; i = next[i]) {
			cap[i] += lambda_dif * cap_b[i];
			if(flow[i] <= cap[i]) continue;
			flow_type val = flow[i] - cap[i];

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
void static_flow_graph::update_lambda(flow_type lambda) {
	flow_type extra_excess = 0;
	for(int i = head[s]; i != -1; i = next[i]) {
    	cap[i] = cap_a[i] + lambda * cap_b[i];
    	if(h[to[i]] >= n) continue;
    	flow_type val = cap[i] - flow[i];
    	extra_excess += val;
    	//ex[to[i]] += val;
    	//ex[s] -= val;
    	flow[i] = cap[i];
    	flow[i^1] = -flow[i];           			
	}

	for(int i = head[t]; i != -1; i = next[i]) {
    	cap[i^1] = cap_a[i^1] + lambda * cap_b[i^1];
    	if(flow[i^1] <= cap[i^1]) continue;
    	flow_type val = flow[i^1] - cap[i^1];
        extra_excess += val;
        //ex[t] -= val;
		//ex[to[i]] += val;
		flow[i^1] = cap[i^1];
		flow[i] = -flow[i^1];
	}
	for(int i = 0; i < n; ++i) {
		ex[i] = 0;
		for(int j = head[i]; j != -1; j = next[j]) {
			ex[i] -= flow[j];						
		}
	}
#ifdef DEBUG
	print("Extra excess = ", 3);
	print(extra_excess, 3);
	print("\n", 3);
#endif
}

void static_flow_graph::push(int v, int e) {
	pushCounter++;
	workDone += PUSH_WORK_CONST;
	if(to[e] != s && to[e] != t && is_zero(ex[to[e]])) {
		rem(inactive[h[to[e]]], to[e]);
		maxheight = max(maxheight, h[to[e]]);
		add_front(active[h[to[e]]], to[e]);
	}
	flow_type val = (ex[v] < cap[e] - flow[e] ? ex[v] : cap[e] - flow[e] );
#ifdef DEBUG
	if(is_zero(val)) {
		print("Warning! Pushing zero value!!!\n", 1);
		print("Current excess of ", 1);
		print(v, 1);
		print(" is ", 1);
		print(ex[v], 1);
		print("\n", 1);
		print("Residual capacity = ", 1);
		print(cap[e] - flow[e], 1);
		print("\nPush was performed from ", 1);
		print(v, 1);
		print(" at height ", 1);
		print(h[v], 1);
		print(" to ", 1);
		print(to[e], 1);
		print("\n", 1);
	}
	print("Current excess of ", 3);
	print(v, 3);
	print(" is ", 3);
	print(ex[v], 3);
	print("\n", 3);
#endif		

#ifdef DEBUG
	print("Pushing ", 3);
	print(val, 3);
	print(" from ", 3);
	print(v, 3);
	print(" to ", 3);
	print(to[e], 3);
	print("\n", 3);
#endif

    ex[to[e]] += val;
	if(ex[v] < cap[e] - flow[e]) {
		flow[e] += val;
		ex[v] = 0;
	} else{
		flow[e] = cap[e];
		ex[v] -= val;
	}
	flow[e^1] = -flow[e];
	if(is_zero(ex[v])) {
		rem(active[h[v]], v);
		add_front(inactive[h[v]], v);
	}	
}

void static_flow_graph::global_relabeling(int stage) {
#ifdef DEBUG	
	print("Starting global relabeling at stage ", 1);
	print(stage + 1, 1);
	print("\nNext update at ", 1);
	print(nextUpdate, 1);
	print("\n", 1);
	//show(3);
#endif
	globalUpdates++;
	queue<int> q;
	if(stage == 0) {
		q.push(t);
	} else {
		q.push(s);
	}
	for(int i = 0; i < n; ++i) {
		g[i] = (i == t) || (i == s);
	}
	for(int i = 0; i < 2 * n; ++i) {
		active[i] = inactive[i] = -1;
	}
	for(int i = 0; i < n; ++i) {
		nextB[i] = -1;
		prevB[i] = -1;
	}	
	maxheight = 0;
	//min_residual = 1e100;
	while(!q.empty()) {
		int temp = q.front();
		q.pop();
		for(int i = head[temp]; i != -1; i = next[i]) {
			if(!g[to[i]] && (!is_zero(cap[i^1] - flow[i^1]))) {
#ifdef DEBUG
				if(h[to[i]] > h[temp] + 1) {
					print("Node ", 1);
					print(to[i], 1);
					print(" has residual arc with cap ", 1);
					print(cap[i^1] - flow[i^1], 1);
					print(" at ", 1);
					print(temp, 1);
					print("\n", 1);
					print(h[to[i]], 1);
					print(" ", 1);
					print(h[temp] + 1, 1);
					print("  BADBADBADB\n", 1);
				}
#endif
				//min_residual = min(min_residual, cap[i^1] - flow[i^1]);

				h[to[i]] = h[temp] + 1;
				g[to[i]] = 1;
				q.push(to[i]);
				if(!is_zero(ex[to[i]])) {
					add_front(active[h[to[i]]], to[i]);
					maxheight = h[to[i]];
				} else {
					add_front(inactive[h[to[i]]], to[i]);
				}
			}
		}
	}
#ifdef DEBUG
	print("Finished global relabeling\n", 1);
	//print_heights(1);
	//print("Minimum capacity of used residual edge = ", 1);
	//print(min_residual, 1);
	//print("\n", 1);
#endif
}


void static_flow_graph::relabel(int v, int stage) {                              	
#ifdef DEBUG	
	print("Relabeling ", 3);
	print(v, 3);
	print("\n", 3);
#endif	
	relabelCounter++;
	workDone += RELABEL_WORK_CONST;
	rem(active[h[v]], v);
	if(active[h[v]] == -1 && inactive[h[v]] == -1) {
#ifdef DEBUG
      	if(h[v] >= n) {
      		print("BADBADBAD!!!! Oh, by the way ex[", 1);
      		print(v, 1);
      		print("]=", 1);
      		print(ex[v], 1);
      		print("\n", 1);
      	}
#endif		
		cur_edge[v] = head[v];
		gap(h[v]);
		h[v] = n + 1;
		return;
	}
	int u = -1;
	for(int i = head[v]; i != -1; i = next[i]) {	
		if(!is_zero(cap[i] - flow[i]) && (u == -1 || h[to[u]] > h[to[i]])) {
			u = i;
		}
		workDone += RELABEL_WORK_PER_ARC;
	}
#ifdef DEBUG
	print("New label is ", 3);
	print(h[to[u]] + 1, 3);
	print("\n", 3);
#endif	
	cur_edge[v] = u;
	h[v] = h[to[u]] + 1;	
	if(stage == 1 || h[v] < n) {
		maxheight = max(maxheight, h[v]);
		add_front(active[h[v]], v);
	}
}

int static_flow_graph::node_for_discharge(int stage) {
	while(maxheight >= stage * n && active[maxheight] == -1) {
		maxheight--;
	}
	if(maxheight == -1) return -1;
#ifdef DEBUG
	if(active[maxheight] != -1 && is_zero(ex[active[maxheight]])) {
		print("Warning!!! Popped for discharge node ", 1);
		print(active[maxheight], 1);
		print(" has excess of ", 1);
		print(ex[active[maxheight]], 1);
		print("\n", 1);
		print_heights(1);
	}
#endif
	return active[maxheight];
}
	
void static_flow_graph::discharge(int v, int stage) {
	if(head[v] == -1) {
#ifdef DEBUG
		print("What?\n", 1);
#endif
		return;
	}
	while(1) {		
#ifdef DEBUG
		print("Discharging ", 2);
		print(v, 2);
		print("\n", 2);
#endif
		if(!is_zero(cap[cur_edge[v]] - flow[cur_edge[v]]) && h[v] == h[to[cur_edge[v]]] + 1) {
			push(v, cur_edge[v]);
		}
		if(is_zero(ex[v])) break;
		cur_edge[v] = next[cur_edge[v]];
		if(cur_edge[v] == -1) {
			relabel(v, stage);
			if(stage == 0 && h[v] >= n) {
				break;
			}
		}
#ifdef DEBUG
		print("Discharge finished\n", 2);
#endif
	}
}

void static_flow_graph::gap(int height) {
#ifdef DEBUG
	print("Starting gap\n", 2);
	//print_heights(3);
#endif	
	gapCounter++;
	for(int i = height + 1; inactive[i] != -1; ++i) {
		/*
		while(active[i] != -1) {
			h[active[i]] = n + 1;
			active[i] = nextB[active[i]];
		}
		*/
		while(inactive[i] != -1) {
			h[inactive[i]] = n + 1;
			rem(inactive[i], inactive[i]);
			//inactive[i] = nextB[inactive[i]];
		}
	}
	//maxheight = height - 1;
#ifdef DEBUG
	print("Finished gap\n", 2);
#endif	
}
	
void static_flow_graph::_add_edge(int a, int b, flow_type A, flow_type B) {
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
   	active = new int[2 * N];
   	inactive = new int[2 * N];
   	for(int i = 0; i < 2 * N; ++i) {
   		active[i] = inactive[i] = -1;
   	}

   	g = new bool[N];
   	cur_edge = new int[N];
   	ex = new flow_type[N];
   	h = new int[N];

    next = new int[2 * M];
   	to = new int[2 * M];
    for(int i = 0; i < N; ++i) {
    	head[i] = -1;
    }
    cap_a = new flow_type[2 * M];
    cap_b = new flow_type[2 * M];
    cap = new flow_type[2 * M];
    flow = new flow_type[2 * M];
        
    reversed = false;
    n = N;
    m = 2 * M;
    s = S;
    t = T;
    edge = 0;    			
	//OBSOLETTE
    eps = 1e-9;
}

void static_flow_graph::add_edge(int a, int b, flow_type A, flow_type B) {
	_add_edge(a, b, A, B);
	_add_edge(b, a, 0, 0);		
}     
    
flow_type static_flow_graph::max_flow(int stage) {
#ifdef DEBUG
	print("Starting maxflow stage ", 3);
	print(stage + 1, 3);
	print("\n", 3);
	//show(1);
#endif	
	while(1) {
		int node = node_for_discharge(stage);
		if(node == -1) break;
#ifdef DEBUG		
		print("Max height = ", 2);
		print(maxheight, 2);
		print("\nDischarge node at height ", 2);
		print(maxheight, 2);
		print(": ", 2);
		print(active[maxheight], 2);
		print("\n", 2);
#endif
		discharge(node, stage);
		if(workDone >= nextUpdate) {
			global_relabeling(stage);
			nextUpdate += globalUpdateBarrier;
		}
	}
#ifdef DEBUG
	print("Maxflow completed!!!\n", 1);
	//show(1);
#endif	
	return ex[t];	
}

flow_type static_flow_graph::max_flow() {
	clear_stats();
	for(int i = 0; i < n; ++i) {
		h[i] = n * (i == s);
		ex[i] = 0;
		cur_edge[i] = head[i];
	}
	for(int i = head[s]; i != -1; i = next[i]) {
		ex[to[i]] = cap[i];
		ex[s] -= cap[i];
		flow[i] = cap[i];
		flow[i^1] = -flow[i];
	}
	maxheight = 0;
#ifdef DEBUG   	
   	print("Stage ", 1);
   	print(0, 1);
   	print("\n", 1);
#endif   	
   	global_relabeling(0);
	nextUpdate += globalUpdateBarrier;
    max_flow(0);

#ifdef DEBUG 
	print("Stage ", 1);
   	print(1, 1);
   	print("\n", 1);
#endif
   	global_relabeling(1);
	nextUpdate += globalUpdateBarrier;
    return max_flow(1);
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
		if(!is_zero(cap[i] - flow[i]) && has_flow(to[i])) return true;
	}
	return false;
}
	

flow_type static_flow_graph::breakpoint_SIMPLE(flow_type init_lambda) {
#ifdef DEBUG
	print("Mode: simple; Stage 1\n", 1);
#endif	
	clear_stats();
	reverse_flow();
	
	initialize_flow(init_lambda);
	for(int i = 0; i < n; ++i) {
		h[i] = n * (i == s);
		ex[i] = 0;
		cur_edge[i] = head[i];
	}
	for(int i = head[s]; i != -1; i = next[i]) {
		ex[to[i]] = cap[i];
		ex[s] -= cap[i];
		flow[i] = cap[i];
		flow[i^1] = -flow[i];
	}
	maxheight = 0;
	lambda = init_lambda;
		
    global_relabeling(0);
	nextUpdate += globalUpdateBarrier;
	flow_type fl = max_flow(0);

	while(1) {
		queue<int> q;
		//reverse_flow();
		for(int i = 0; i < n; i++) {
			g[i] = (i == t);
		}			
		q.push(t);
		
		int cut_size = 0;
		//flow_type min_residual = 1000000000;

		while(!q.empty()) {
			cut_size++;
			int temp = q.front();
			//cerr << "q:" << temp << endl;
			q.pop();
			for(int i = head[temp]; i != -1; i = next[i]) {
				if(!is_zero(cap[i^1] - flow[i^1]) && !g[to[i]]) {
					//min_residual = min(min_residual, cap[i^1] - flow[i^1]);
					q.push(to[i]);
					g[to[i]] = true;
				}
			}
		}
#ifdef DEBUG
		//print("Minimum capacity of used residual edge = ", 1);
		//print(min_residual, 1);
		//print("\n", 1);
		print("Cut size == ", 1);
		print(cut_size, 1);
		print("\n", 1);
#endif
		if(cut_size == 1) {
			//reverse_flow();
			break;
		}
	    g[t] = false;
		flow_type mult = 0, stat = 0;
		for(int i = 0; i < n; ++i) {
			if(g[i]) {
				for(int j = head[i]; j != -1; j = next[j]) {
					if(!g[to[j]]) {
						stat -= cap_a[j^1];
						mult += cap_b[j^1];
					}
				}
			}
			if(i == t) {
				for(int j = head[i]; j != -1; j = next[j]) {
					if(g[to[j]]) {
						stat += cap_a[j^1];
						mult -= cap_b[j^1];
					}
				}
			}
		}
		//reverse_flow();

#ifdef DEBUG
		print(stat, 1);
		print(" ", 1);
		print(mult, 1);
		print("\n", 1);
#endif
		flow_type new_lambda = stat / mult;
#ifdef DEBUG
		print("new lambda = ", 1);
		print(new_lambda, 1);
		print("\n", 1);
#endif		
		//cerr << s_load << " " << cut << " " << tau << " " << new_tau << endl;
	    //decrease_cap(new_lambda - lambda);
	    if(is_zero(new_lambda - lambda)) {
	    	break;
	    }
	    //update_lambda(new_lambda);    
	    
	    lambda = new_lambda;
	    initialize_flow(lambda);
        for(int i = 0; i < n; ++i) {
			h[i] = n * (i == s);
			ex[i] = 0;
			cur_edge[i] = head[i];
		}
		for(int i = head[s]; i != -1; i = next[i]) {
			ex[to[i]] = cap[i];
			ex[s] -= cap[i];
			flow[i] = cap[i];
			flow[i^1] = -flow[i];
		}	

	   	global_relabeling(0);
		nextUpdate += globalUpdateBarrier;
		fl = max_flow(0);
		//cerr << 1 / lambda << " " << fl << endl;
	}
#ifdef DEBUG
	print("Stage 2\n", 1);
#endif	
	global_relabeling(1);
	nextUpdate += globalUpdateBarrier;
	max_flow(1);
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
flow_type static_flow_graph::breakpoint_GGT(flow_type init_lambda) {
#ifdef DEBUG
	print("Mode: GGT; Stage 1\n", 1);
#endif	
	initialize_flow(init_lambda);
	reverse_flow();
	
	for(int i = 0; i < n; ++i) {
		h[i] = n * (i == s);
		ex[i] = 0;
		cur_edge[i] = head[i];
	}
	for(int i = head[s]; i != -1; i = next[i]) {
		ex[to[i]] = cap[i];
		ex[s] -= cap[i];
		flow[i] = cap[i];
		flow[i^1] = -flow[i];
	}
	maxheight = 0;
	lambda = init_lambda;
		
    global_relabeling(0);
	nextUpdate += globalUpdateBarrier;
	flow_type fl = max_flow(0);

	while(1) {
		queue<int> q;
		//reverse_flow();
		for(int i = 0; i < n; i++) {
			g[i] = (i == t);
		}			
		q.push(t);
		
		int cut_size = 0;
		//flow_type min_residual = 1e100;

		while(!q.empty()) {
			cut_size++;
			int temp = q.front();
			//cerr << "q:" << temp << endl;
			q.pop();
			for(int i = head[temp]; i != -1; i = next[i]) {
				if(!is_zero(cap[i^1] - flow[i^1]) && !g[to[i]]) {
					//min_residual = min(min_residual, cap[i^1] - flow[i^1]);
					q.push(to[i]);
					g[to[i]] = true;
				}
			}
		}
#ifdef DEBUG
		//print("Minimum capacity of used residual edge = ", 1);
		//print(min_residual, 1);
		//print("\n", 1);
		print("Cut size == ", 1);
		print(cut_size, 1);
		print("\n", 1);
#endif
		if(cut_size == 1) {
			//reverse_flow();
			break;
		}
	    g[t] = false;
		flow_type mult = 0, stat = 0;
		for(int i = 0; i < n; ++i) {
			if(g[i]) {
				for(int j = head[i]; j != -1; j = next[j]) {
					if(!g[to[j]]) {
						stat -= cap_a[j^1];
						mult += cap_b[j^1];
					}
				}
			}
			if(i == t) {
				for(int j = head[i]; j != -1; j = next[j]) {
					if(g[to[j]]) {
						stat += cap_a[j^1];
						mult -= cap_b[j^1];
					}
				}
			}
		}
		//reverse_flow();

#ifdef DEBUG
		print(stat, 1);
		print(" ", 1);
		print(mult, 1);
		print("\n", 1);
#endif
		flow_type new_lambda = stat / mult;
#ifdef DEBUG
		print("new lambda = ", 1);
		print(new_lambda, 1);
		print("\n", 1);
#endif		
		//cerr << s_load << " " << cut << " " << tau << " " << new_tau << endl;
	    //decrease_cap(new_lambda - lambda);
	    if(is_zero(new_lambda - lambda)) {
	    	break;
	    }
	    update_lambda(new_lambda);    
	    lambda = new_lambda;
	    
	   	global_relabeling(0);
		nextUpdate += globalUpdateBarrier;
		fl = max_flow(0);
		//cerr << 1 / lambda << " " << fl << endl;
	}
#ifdef DEBUG
	print("Stage 2\n", 1);
#endif	
	global_relabeling(1);
	nextUpdate += globalUpdateBarrier;
	max_flow(1);
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

flow_type static_flow_graph::leftmost_breakpoint(flow_type init_lambda) {
#ifdef DOUBLE_FLOW	
	return breakpoint_GGT(init_lambda);
#endif
#ifdef RATIONAL_FLOW
	return breakpoint_SIMPLE(init_lambda);
#endif
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
		print(", height = ", priority);
		print(h[i], priority);
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
		if(temp == -1) continue;
		print("Active vertices at height ", priority);
		print(i, priority);
		print(":", priority);
		while(temp != -1) {
			print(temp, priority);
			print(" (ex=", priority);
			print(ex[temp], priority);
			print(") ", priority);
			temp = nextB[temp];
		}
		print(" (0-indexed)\n", priority);
	}
	for(int i = 0; i < 2 * n; ++i) {
		int temp = inactive[i];
		if(temp == -1) continue;
		print("Inactive vertices at height ", priority);
		print(i, priority);
		print(":", priority);
		while(temp != -1) {
			print(temp, priority);
			print(" (ex=", priority);
			print(ex[temp], priority);
			print(") ", priority);
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
	print("Nubmer of gaps: ", priority);
	print(gapCounter, priority);
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

flow_type static_flow_graph::get_flow(int idx) {
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

long long static_flow_graph::gaps() {
	return gapCounter;
}
