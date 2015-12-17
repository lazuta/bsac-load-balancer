#include <queue>
#include <cstdlib>
#include <cmath>

using namespace std;


class static_flow_graph {
private:	
	int *head, *next, *to, *prev, *h, *cur_edge;
	double *cap, *cap_a, *cap_b, *flow, *ex, *lambda_i;
	bool *g;
	bool reversed;
	int edge, n, m, s, t;
	double eps;
	double lambda;
	queue<int> active_vrtx;

	void initialize_flow(double L) {
		lambda = L;
		for(int i = 0; i < m; ++i) {
			cap[i] = cap_a[i] + lambda * cap_b[i];
			flow[i] = 0;
		}
	}

	/**
	 * lambda_dif should be non-positive.
	 */
	void decrease_cap(double lambda_dif) {
		for(int v = 0; v < n; ++v) {
			for(int i = head[v]; i != -1; i = next[i]) {
				cap[i] += lambda_dif * cap_b[i];
				if(flow[i] <= cap[i]) continue;
				double val = flow[i] - cap[i];
				if(to[i] == t || v == s) {
                	if(ex[v] < eps) {
						active_vrtx.push(v);
					}
				}

                ex[v] += val;
				ex[to[i]] -= val;
				flow[i] -= val;
				flow[i^1] = -flow[i];
			}
		}
	}

	void push(int v, int e) {
		if(to[e] != s && to[e] != t && ex[to[e]] == 0) {
			active_vrtx.push(to[e]);
		}
		double val = min(ex[v], cap[e] - flow[e]);
        ex[v] -= val;
		ex[to[e]] += val;
		flow[e] += val;
		flow[e^1] = -flow[e];	
	}

	void relabel(int v) {                              	
		int u = -1;
		for(int i = head[v]; i != -1; i = next[i]) {	
			if(cap[i] - flow[i] > 0 && (u == -1 || h[u] > h[to[i]])) {
				u = to[i];
			}
		}
		h[v] = h[u] + 1;	
	}
	
	void discharge(int v) {
		if(head[v] == -1) return;
		while(1) {
			if(cap[cur_edge[v]] > flow[cur_edge[v]] && h[v] == h[to[cur_edge[v]]] + 1) {
				push(v, cur_edge[v]);
			}
			if(ex[v] == 0) break;
			
			cur_edge[v] = next[cur_edge[v]];
			if(cur_edge[v] == -1) {
				relabel(v);
				cur_edge[v] = head[v];
			}
		}
	}
	
	void _add_edge(int a, int b, double A, double B) {
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


public:
	void initialize_graph(int N, int M, int S, int T, double e) {
    	head = new int[N];
    	prev = new int[N];
    	g = new bool[N];
    	cur_edge = new int[N];
    	ex = new double[N];
    	h = new int[N];
    	lambda_i = new double[N];

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

        eps = e;
	}

	void add_edge(int a, int b, double A, double B) {
		_add_edge(a, b, A, B);
		_add_edge(b, a, 0, 0);		
	}     
    
    double max_flow() {
        while(!active_vrtx.empty()) {
			discharge(active_vrtx.front());
			active_vrtx.pop();
		}
		return ex[t];	
	}

	
	/**
	 * Inverses current flow state. That is, edges and their capacity and flow are reoriented 
	 * backwards, sink and source nodes are swapped.
	 */
	void reverse_flow() {
		reversed = !reversed;
		swap(s, t);
		for(int i = 0; i < m / 2; ++i) {
			swap(cap[2 * i], cap[2 * i + 1]);
			swap(cap_a[2 * i], cap_a[2 * i + 1]);
            swap(cap_b[2 * i], cap_b[2 * i + 1]);
            swap(flow[2 * i], flow[2 * i + 1]);
		} 
	}
	
	bool has_flow(int v) {
		if(v == t) return true;
		if(g[v]) return false;
		g[v] = true;

		for(int i = head[v]; i != -1; i = next[i]) {
			if(cap[i] - flow[i] > eps && has_flow(to[i])) return true;
		}
		return false;
	}
	
	/**
	 *
	 */
	double leftmost_breakpoint_2(double init_lambda) {
		initialize_flow(init_lambda);
		double load = 0;
		for(int i = head[s]; i != -1; i = next[i]) {
        	load += cap[i];
        }

        cerr << "Init labmda " << init_lambda << endl;
		reverse_flow();
		
		for(int i = 0; i < n; ++i) {
			h[i] = n * (i == s);
			ex[i] = 0;
			cur_edge[i] = head[i];
		}
		for(int i = head[s]; i != -1; i = next[i]) {
			active_vrtx.push(to[i]);
			ex[to[i]] = cap[i];
			flow[i] = cap[i];
			flow[i^1] = -flow[i];
		}
		for(int i = 0; i < n; ++i) {
			lambda_i[i] = init_lambda;
		}
			    
        
        bool good;
        double l1 = init_lambda, l2 = 0;
        double fl = 0;
        while(1) {
			l1 = init_lambda, l2 = 0;
        	fl = max_flow();
        	reverse_flow();
        	cerr << "Rev = " << reversed << endl;
        	for(int i = 0; i < n; ++i ) {
        		if(i != s && i != t) {
        			cerr << "Lambda[" << i << "] = " << lambda_i[i] << endl;
        		}
        	}
        	for(int i = 0; i < n; ++i) {
        		if(i == s) {
        			for(int j = head[i]; j != -1; j = next[j]) {
        				cerr << "lambda * q[" << to[j] << "] = " << cap[j] << ", flow = " << flow[j] << endl;
        			}
        		}
        		else if(i != t) {
                	for(int j = head[i]; j != -1; j = next[j]) {
                		if(to[j] != t) cerr << "cap[" << i << ", " << to[j] << "] = " << cap[j] << ", flow = " << flow[j] << endl;
                		else cerr << "p[" << i << "] = " << cap[j] << ", flow = " << flow[j] << endl;
                	}        			
        		}
        	}
        	
        	cerr << "Load and flow " << load << " " << fl << endl;
			good = 0;
        	for(int i = head[s]; i != -1; i = next[i]) {
        		if(cap[i] - flow[i] > eps) {
        			cerr << to[i] << " is picked: " << cap[i] - flow[i] << endl;
        			double mult = 0, stat = 0;
					for(int j = head[to[i]]; j != -1; j = next[j]) {
						if((to[j] != s) && (cap[j] - flow[j] <= 0)) {
							stat -= cap_a[j];
							mult += cap_b[j];
						}
						if(to[j] == s) {
							stat += cap_a[j^1];
							mult -= cap_b[j^1];
						}
					
					}
					cerr << to[i] << " " << stat << " " << mult << " " << lambda_i[to[i]] << " " << stat / mult << endl;			
					double new_lambda = stat / mult;
					///------------------
					for(int j = head[to[i]]; j != -1; j = next[j]) {
						if(to[j] == s) {
							cap[j] = cap_a[j] + new_lambda * cap_b[j];
							cap[j^1] = cap_a[j^1] + new_lambda * cap_b[j^1];
							if(flow[j^1] - cap[j^1] < eps) break;
							double val = flow[j^1] - cap[j^1];
						    ex[to[i]] += val;
							ex[s] -= val;
							flow[j^1] -= val;
							flow[j] = -flow[j^1];
							break;
			    		}
					}
					active_vrtx.push(to[i]);
					///------------------
					lambda_i[to[i]] = stat / mult;
					cerr  << "New lambda for " << to[i] << " is " << lambda_i[to[i]] << endl;
					for(int j = head[to[i]]; j != -1; j = next[j]) {
        				if(to[j] == s) cerr << "lambda * q[" << to[i] << "] = " << cap[j^1] << ", flow = " << flow[j^1] << endl;
        			}
        		
					reverse_flow();
					good = true;
					break; 		
        		}
        	}
        	if(!good) break;
        	for(int i = 0; i < n; ++i) {
				if(i != s && i != t) {
					l1 = min(l1, lambda_i[i]);
					l2 = max(l2, lambda_i[i]);
         	   }
			}
			cerr << l1 << " " << l2 << endl;	
		
        }
        
        l1 = init_lambda, l2 = 0;
        for(int i = 0; i < n; ++i) {
			if(i != s && i != t) {
				l1 = min(l1, lambda_i[i]);
				l2 = max(l2, lambda_i[i]);
            }
		}	
		cerr << "Min lambda = " << l1 << endl;
		cerr << "Max lambda = " << l2 << endl;

		if(reversed) reverse_flow();
		return l1;
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
	double leftmost_breakpoint(double init_lambda) {
		initialize_flow(init_lambda);
		reverse_flow();
		
		for(int i = 0; i < n; ++i) {
			h[i] = n * (i == s);
			ex[i] = 0;
			cur_edge[i] = head[i];
		}
		for(int i = head[s]; i != -1; i = next[i]) {
			active_vrtx.push(to[i]);
			ex[to[i]] = cap[i];
			flow[i] = cap[i];
			flow[i^1] = -flow[i];
		}
		lambda = init_lambda;
		
		
        double fl = max_flow();
        double load = 0;
        for(int i = head[s]; i != -1; i = next[i]) {
        	load += cap[i];
        }
		cerr << 1 / lambda << " " << fl << " " << load << endl;
		
		while(fl - load < -eps) {
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
			
			double new_lambda = stat / mult;
    	    //cerr << s_load << " " << cut << " " << tau << " " << new_tau << endl;
    	    //initialize_flow(n, m, new_tau);
    	    decrease_cap(new_lambda - lambda);
    	    lambda = new_lambda;
    	    
    	    fl = max_flow();
    	    cerr << 1 / lambda << " " << fl << endl;
		}
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
	

	double get_flow(int idx) {
		return flow[2 * idx];
	}
};