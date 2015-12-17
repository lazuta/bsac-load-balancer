#include <queue>
#include <cstdlib>
#include <cmath>
//#include <iostream>

using namespace std;

class static_flow_graph {
private:	
	int *head, *next, *to, *prev;
	double *cap, *cap_a, *cap_b, *flow;
	bool *g;
	int edge, n, m, s, t;
	double eps;
	double lambda;

	double add_flow_bfs() {
		queue<int> q;
		for(int i = 0; i < n; i++) {
			prev[i] = -1;
		}
		q.push(s);
		while(!q.empty()) {
			int temp = q.front();
			q.pop();
			for(int i = head[temp]; i != -1; i = next[i]) {
				if(cap[i] - flow[i] > eps && prev[to[i]] == -1) {
					q.push(to[i]);
					prev[to[i]] = i;
				}
			}
		}
		if(prev[t] == -1) return 0;
		double min_cap = 1e100;
		int temp = t;
		while(temp != s) {
			min_cap = min(min_cap, cap[prev[temp]] - flow[prev[temp]]);
			temp = to[prev[temp]^1];
		}
		temp = t;
		while(temp != s) {
			flow[prev[temp]] += min_cap;
			flow[prev[temp]^1] -= min_cap;
			temp = to[prev[temp]^1];
		}
		return min_cap;
	}

public:
	void initialize_graph(int N, int M, int S, int T, double e) {
    	head = new int[N];
    	prev = new int[N];
    	g = new bool[N];

    	next = new int[2 * M];
    	to = new int[2 * M];
    	for(int i = 0; i < N; ++i) {
    		head[i] = -1;
    	}
    	cap_a = new double[2 * M];
    	cap_b = new double[2 * M];
    	cap = new double[2 * M];
        flow = new double[2 * M];
        
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

	void initialize_flow(double L) {
		lambda = L;
		for(int i = 0; i < m; ++i) {
			cap[i] = cap_a[i] + lambda * cap_b[i];
			flow[i] = 0;
		}
	}
	
	void update_cap(double lambda_dif) {
		for(int i = 0; i < m; ++i) {
			cap[i] += lambda_dif * cap_b[i];
		}
	}

	double max_flow() {
		double f = add_flow_bfs();
		double sum = 0;
	    while(f > 0) {
	    	sum += f;
	       	f = add_flow_bfs(); 
	    }
	    return sum;    
	}

	/**
	 * Note that it actually finds the rightmost breakpoint for
	 * particular case when all capacity multipliers have positive value
	 * and the rightmost (with lambda -> +infinity) minimum cut is
	 * {s}|V\{s}
	 *
	 * Positivity of multipliers is a necessary condition for a Newton's like method
	 * used in this function to generate an increasing sequence of parameters that
	 * converge to the rightmost breakpoint.
	 *
	 * Also, init_lambad is required to be positive non-zero and less then rightmost 
	 * breakpoint. Warning!!! Function might not work as expected with init_lambda 
	 * values close to zero (the closer the value to zero the more impact calculation 
	 * errors will have).
	 */
	double rightmost_breakpoint(double init_lambda) {
		initialize_flow(init_lambda);
		lambda = init_lambda;
        		
		double fl = max_flow();

        double load = 0;
        for(int i = head[s]; i != -1; i = next[i]) {
        	load += cap[i];
        }
		//cerr << lambda << " " << fl << " " << load << endl;
		
		while(fl - load < -eps) {
			//cerr << "lambda" << lambda << endl;
			queue<int> q;
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
			/*
			rep(i, 2 * n + m) {
				cerr << cap[2 * i] << " " << flow[2 * i] << endl;
			}
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
			
			double new_lambda = stat / mult;
    	    //cerr << new_lambda  << " " << lambda << "!!!" << endl;
    	    //initialize_flow(n, m, new_tau);
    	    update_cap(new_lambda - lambda);
    	    lambda = new_lambda;
    	    
    	    fl = max_flow();
    	    //cerr << lambda << " " << fl << endl;
		}
		return lambda;
	}
	
	double get_flow(int idx) {
		return flow[2 * idx];
	}
};