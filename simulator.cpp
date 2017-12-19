#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <list>
#include <cstring>

using namespace std;

#include "simulator.h"
#include "utils/logger.h"
#include "utils/rational.h"


bool scheduler::is_saturated() {
   	for(int i = 0; i < saturation_ratio.size(); ++i) {
   		if(saturation_ratio[i].second->scheduled_load < saturation_ratio[i].second->transferred_load ||
   			saturation_ratio[i].second->current_task != NULL) {
   			return false;
 		}
	}
	return true;
}

scheduler::scheduler() {
	task_bound = -1;
	excess = 0;
}

task::task() {
	id = 0;
	processing_time = task_processing_time_expectation + task_processing_time_variance - (rand() % (2 * task_processing_time_variance + 1));
	content_size = task_content_size_expectation + task_content_size_variance - (rand() % (2 * task_content_size_variance + 1));
	current_pu = NULL;
	resolved = false;
}

channel::channel() {
	id = 0;
	in = NULL;
	out = NULL;
	current_progress = 0;
	capacity_ptu = 0;
	delay = 0;
	current_task = NULL;
	scheduled_load = 0;
	transferred_load = 0;
}

list<task*> delayer::get_current() {
	return time_slots[current_time];
}

void delayer::tick() {
	cnt -= time_slots[current_time].size();
	time_slots[current_time].clear();
	current_time = (current_time + 1) % max_delay;
}

void delayer::add_task(task* t, int delay) {
	cnt++;
	time_slots[(current_time + delay) % max_delay].push_back(t);
}

delayer::delayer() {
   	current_time = 0;
   	cnt = 0;
}

bool delayer::empty() {
   	return cnt == 0;    	
}   

int delayer::size() {
   	return cnt;
}                      

processing_unit::processing_unit() {
	id = 0;
	perfomance_ptu = 0;
	current_progress = 0;
	current_task = NULL;
	sch = NULL;
	abstract_delayer = NULL;
}

void processing_unit::add_task(task* tsk) {
	buffer.push(tsk);
}

void processing_unit::add_channel(channel* ch) {
	outgoing_channels.push_back(ch);
	sch->saturation_ratio.push_back(make_pair(1, ch));
}

void read_graph(static_flow_graph* graph, FILE* f) {
	int s, t;
	fscanf(f, "%d%d", &n, &m);
	s = n, t = n + 1;
	task_cnt = 0;
	int num = 0;
	sum_perf = 0;
	graph->initialize_graph(n + 2, 2 * n + m, s, t);
	
	m_t = 0;

	for(int i = 0; i < n; ++i) {
		int p, q;
		fscanf(f, "%d%d", &p, &q);
		processing_unit* u = new processing_unit();
		u->sch = new scheduler();
		u->id = i;
		u->abstract_delayer = new delayer();
		units.push_back(u);
		u->perfomance_ptu = p * task_processing_time_expectation;
		sum_perf += p;
		double load = 0;
		print("Unit ", 4);
		print(i, 4);
		print(" is assigned with ", 4);
		print(q, 4);
		print(" tasks, id's: ", 4);
		print(task_cnt, 4);
		print("-", 4);
		print(task_cnt + q - 1, 4);
		print("\n", 4);
		task_cnt += q;
		//if(q > 0) {
			graph->add_edge(s, i, 0, q);
		//}
		if(simulation != 0) {
        	for(int j = 0; j < q; ++j) {
				task* tsk = new task();
				tsk->current_pu = u;
				tsk->id = num++;
				u->add_task(tsk);
				load += tsk->processing_time;
			}
		}
		m_t = m_t >  load / u->perfomance_ptu ? m_t : load / u->perfomance_ptu;
	}
	for(int i = 0; i < m; ++i) {
		int q, w; double c;
		fscanf(f, "%d%d%lf", &q, &w, &c);
		
		channel* ch = new channel();
		ch->in = units[q - 1];
		ch->out = units[w - 1];
		ch->id = i;
		ch->delay = 1;
		ch->capacity_ptu = c * task_content_size_expectation;

		units[q - 1]->add_channel(ch);
	}	

	for(int i = 0; i < n; ++i) {
		graph->add_edge(i, t, units[i]->perfomance_ptu / task_processing_time_expectation, 0);
    }

    for(int i = 0; i < n; ++i) {
    	for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
		    graph->add_edge(i, units[i]->outgoing_channels[j]->out->id, units[i]->outgoing_channels[j]->capacity_ptu / task_content_size_expectation, 0);
    	}
    }

    for(int i = 0; i < n; ++i) {
    	//graph->add_edge(s, i, 0, units[i]->buffer.size());
    }		
}

void calculate_flow(static_flow_graph* graph) {
    print("Flow calculation ...\n", 2);
    scheduling_time = clock();
    //sum_perf = 4;
    //task_cnt = 993;
#ifdef DOUBLE_FLOW
    tau = 1 / graph->leftmost_breakpoint((double)sum_perf /  task_cnt);
#endif
#ifdef RATIONAL_FLOW
	tau = 1 / graph->leftmost_breakpoint(Rational(sum_perf, task_cnt));
#endif
    print("Done\n", 2);
	print("Flow calculation took ", 2);
	print((double)(scheduling_time = clock() - scheduling_time) / CLOCKS_PER_SEC, 2);
	print("sec.\n", 2);
	graph->print_stats(2);
	
	print("Found tau = ", 2);
	print(tau, 2);
	print("\n", 2);

#ifdef DEBUG
	graph->show(3);
#endif
}

void tick(processing_unit* unit) {
	int prev = unit->current_progress;
	unit->current_progress += unit->perfomance_ptu;
	if(unit->current_task != NULL && unit->current_progress >= unit->current_task->processing_time) {
		unit->current_progress -= unit->current_task->processing_time;
		unit->current_task->resolved = true;
		task_resolved++;
		print("Task #", 4);
		print(unit->current_task->id, 4);
		print(" is resolved at ", 4);
		print(unit->id, 4);
		print(" (operations done: ", 4);
		print(unit->current_task->processing_time - prev, 4);
		print(").\n", 4);
		//cerr << "Task #" << unit->current_task->id << " is resolved!!! Congratulations!!!" << endl;
		unit->current_task = NULL;
	}
	if(unit->current_task == NULL) {
		while(!unit->buffer.empty() && unit->current_progress >= unit->buffer.front()->processing_time) {
			unit->current_progress -= unit->buffer.front()->processing_time;
			unit->buffer.front()->resolved = true;
			task_resolved++;
			print("Task #", 4);
			print(unit->buffer.front()->id, 4);
			print(" is resolved at ", 4);
			print(unit->id, 4);
			print(" (operations done: ", 4);
			print(unit->buffer.front()->processing_time, 4);
			print(").\n", 4);
			//cerr << "Task #" << unit->buffer.front()->id << " is resolved!!! Congratulations!!!" << endl;
			unit->buffer.pop();    				
		}
		if(unit->buffer.empty()) {
			unit->current_progress = 0;
		} else {
			unit->current_task = unit->buffer.front();
			print("Task #", 4);
			print(unit->current_task->id, 4);
			print(" is currently being resolved at ", 4);
			print(unit->id, 4);
			print(" (current progress: ", 4);
			print(unit->current_progress, 4);
			print(").\n", 4);
			unit->buffer.pop();
		}
	}                   	
	
	
	if(algorithm == CONSENSUS) {
		clock_t clck = clock();
		for(int i = 0; i < unit->sch->saturation_ratio.size(); ++i) {
    		channel* ch = unit->sch->saturation_ratio[i].second;
    		ch->transferred_load = 0;
    		
			if(unit->buffer.size() * ch->out->perfomance_ptu > 
    			unit->perfomance_ptu * ch->out->buffer.size()) {
    			ch->scheduled_load = alpha * (- ch->out->buffer.size() / ch->out->perfomance_ptu
    				+ unit->buffer.size() / unit->perfomance_ptu);
    			unit->sch->saturation_ratio[i].first = (double)ch->transferred_load / ch->scheduled_load;        		
			} else {
    			ch->scheduled_load = 0;
    			unit->sch->saturation_ratio[i].first = 1;
			}
		//if(ch->scheduled_load) cerr << ch->scheduled_load << "is scheduled to transmit from " << unit->id  << " to " << ch->out->id << endl;
        }
    	scheduling_time += (clock() - clck);            	
		sort(unit->sch->saturation_ratio.begin(), unit->sch->saturation_ratio.end());		
	}

	if(algorithm == ADAPTIVE_FLOW) {
		clock_t clck = clock();
		for(int i = 0; i < unit->sch->saturation_ratio.size(); ++i) {
    		channel* ch = unit->sch->saturation_ratio[i].second;
    		ch->transferred_load = 0;
    		
			long long new_f;
			if(ch->out->buffer.size() + unit->buffer.size() == 0) {
				new_f = 0;
			} else {
				new_f = (-ch->out->sch->excess * unit->buffer.size() * unit->buffer.size() + ch->out->buffer.size() * ch->out->buffer.size() * unit->sch->excess) 
					/ ((long long)ch->out->buffer.size() * (long long)ch->out->buffer.size() + (long long)unit->buffer.size() * (long long)unit->buffer.size());
				if(ch->scheduled_load + new_f * task_content_size_expectation < 0) new_f = -ch->scheduled_load / task_content_size_expectation;
				if(ch->scheduled_load + new_f * task_content_size_expectation > ch->capacity_ptu) 
					new_f = (ch->capacity_ptu - ch->scheduled_load) / task_content_size_expectation;
			}
			ch->out->sch->excess += /*ch->scheduled_load / task_content_size_expectation*/ new_f;
			unit->sch->excess -= /*ch->scheduled_load / task_content_size_expectation*/ new_f;
			ch->scheduled_load += new_f * task_content_size_expectation;

		//if(ch->scheduled_load) cerr << ch->scheduled_load << "is scheduled to transmit from " << unit->id  << " to " << ch->out->id << endl;
        }
    	scheduling_time += (clock() - clck);            	
		sort(unit->sch->saturation_ratio.begin(), unit->sch->saturation_ratio.end());
	}
	
	for(int i = 0; i < unit->sch->saturation_ratio.size(); ++i) {
		channel* ch = unit->sch->saturation_ratio[i].second;
		prev = ch->current_progress;
		ch->current_progress += ch->capacity_ptu;
		if(ch->current_task != NULL && ch->current_progress >= ch->current_task->content_size) {
			/*if(ch->scheduled_load - ch->transferred_load < ch->current_task->content_size - ch->current_progress) {
				ch->scheduled_load = ch->current_task->content_size - ch->current_progress + ch->transferred_load;
			}*/
			ch->current_progress -= ch->current_task->content_size;
			print("Task #", 4);
			print(ch->current_task->id, 4);
			print(" is sent out from ", 4);
			print(ch->in->id, 4);
			print(" (bits sent: ", 4);
            print(ch->current_task->content_size, 4);
			print(").\n", 4);
			ch->out->abstract_delayer->add_task(ch->current_task, ch->delay);
			ch->current_task = NULL;
		}                                       
		if(ch->current_task == NULL) {
			while(!unit->buffer.empty() && ch->transferred_load < ch->scheduled_load && ch->current_progress >= unit->buffer.front()->content_size) {
				ch->current_progress -= unit->buffer.front()->content_size;
				ch->transferred_load += unit->buffer.front()->content_size;
				print("Task #", 4);
				print(unit->buffer.front()->id, 4);
				print(" is sent out from ", 4);
				print(unit->id, 4);
				print(" at ", 4);
				print(ch->out->id, 4);
				print(" (bits sent: ", 4);
                print(unit->buffer.front()->content_size, 4);
				print(").\n", 4);
			    ch->out->abstract_delayer->add_task(unit->buffer.front(), ch->delay);
				unit->buffer.pop();
			}
			if(unit->buffer.empty() || ch->transferred_load >= ch->scheduled_load) {
				ch->current_progress = 0;
			} else {
				ch->current_task = unit->buffer.front();
				print("Task #", 4);
				print(ch->current_task->id, 4);
				print(" is currently being transferred out from ", 4);
				print(unit->id, 4);
				print(" at ", 4);
				print(ch->out->id, 4);
				print(" (current progress: ", 4);
				print(ch->current_progress, 4);
				print(").\n", 4);
				ch->transferred_load += unit->buffer.front()->content_size;
				unit->buffer.pop();
			}
			if(ch->scheduled_load == 0) {
				unit->sch->saturation_ratio[i].first = 100;
			} else {
				unit->sch->saturation_ratio[i].first = (double)ch->transferred_load / ch->scheduled_load;        		
			}
		}
	}
	sort(unit->sch->saturation_ratio.begin(), unit->sch->saturation_ratio.end());	
	list<task*> tmp = unit->abstract_delayer->get_current();

	while(!tmp.empty()) {
		unit->buffer.push(tmp.front());
		print("Task #", 4);
		print(tmp.front()->id, 4);
		print(" is recieved on ", 4);
		print(unit->id, 4);
		print(" (bits recieved: ", 4);
		print(tmp.front()->content_size, 4);
		print(").\n", 4);
		//cerr << "Task #" << tmp.front()->id << " is transferred to unit #" << unit->id << endl;
		tmp.pop_front();
	}
	unit->abstract_delayer->tick();

	unit->status = ((int)(!unit->buffer.empty() || (unit->current_task != NULL) || (!unit->abstract_delayer->empty())) << 1) + (!unit->sch->is_saturated());
	//cerr << "Unit #" << unit->id << " status is " << unit->status << endl;
}

bool tick() {
	bool done = true;
	if(algorithm == ADAPTIVE_FLOW) {
	for(int i = 0; i < units.size(); ++i) {
		print("Excess of node #", 4);
		print(i, 4);
		print(" = ", 4);
		print(units[i]->sch->excess, 4);
		print(" // ", 4);
		print(units[i]->buffer.size(), 4);
		print("\n", 4);
	}		
	}
	for(int i = 0; i < units.size(); ++i) {
		tick(units[i]);
		if(units[i]->status > 1) done = false;
	}

	return done;	
}   

void set_consensus_step(double step) {
	alpha = step;
}

void set_balancing_algorithm(char* algorithm_name) {
	if(strcmp(algorithm_name, "adaptive_flow") == 0) {
		algorithm = ADAPTIVE_FLOW;
	} else if(strcmp(algorithm_name, "consensus") == 0) {
		algorithm = CONSENSUS;
	} else {
		algorithm = PARAMETRIC_FLOW;
	}
}

void set_time_output_step(int step) {
	if(step <= 0) time_output_step = 1;
	else time_output_step = step;
}

void set_simulation(int Do) {
	simulation = Do;
}

void set_processing_time_e(int t) {
	task_processing_time_expectation = t;
}

void set_processing_time_v(int t) {
	task_processing_time_variance = t;
}

void set_content_size_e(int t) {
	task_content_size_expectation = t;
}

void set_content_size_v(int t) {
	task_content_size_variance = t;
}	


void simulate(string path) {
	FILE* _f;
	print(path, 4);
	print("\n", 4);	
	if( ! (_f = (fopen(path.c_str(), "r"))) ) {
		//cout << "Test doesn't exist" << endl;
		print("Test doesn't exist\n", 0);		
		return;
	}	

	print("Balancing algorithm in use: ", 2);
	if(algorithm == PARAMETRIC_FLOW) {
		print("parametric flow\n", 2);
	} else if (algorithm == CONSENSUS) {
		print("consensus\n", 2);
		print("Step size is ", 2);
		print(alpha, 2);
		print("\n", 2);
	} else {
		print("adaptive flow\n", 2);
	}                  	

	static_flow_graph* graph;
	graph = new static_flow_graph();
    read_graph(graph, _f);

	print("Milestone 1\n", 4);

	print("Number of processing units: ", 4);
	print(static_cast<int>(units.size()), 4);
	print("\n", 4);

	print("Number of communication channels: ", 4);
	print(m, 4);
	print("\n", 4);

    //graph->show(4);
    if(algorithm == PARAMETRIC_FLOW) {
    	calculate_flow(graph);
	
		int num = 0;
    	for(int i = 0; i < n; ++i) {
    		for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
    			print("Setting flow on (", 4);
    			print(units[i]->id, 4);
    			print(", ", 4);
    			print(units[i]->outgoing_channels[j]->id, 4);
				print(") = ", 4);
				double f; 
#ifdef RATIONAL_FLOW
				f = (graph->get_flow(2 * n + num) * tau).toDouble();
#endif
#ifdef DOUBLE_FLOW
				f = (graph->get_flow(2 * n + num) * tau);
#endif                                          				
				print(ceil(f) * task_content_size_expectation, 4);
				print("\n", 4);
    			units[i]->outgoing_channels[j]->scheduled_load = ceil(f) * task_content_size_expectation;
    			//cerr << "Scheduled load: " << units[i]->outgoing_channels[j]->scheduled_load << " " << graph->get_flow(n + num) * tau << endl;
    			num++; 
    		}
    	}       
    }	
    if(algorithm == ADAPTIVE_FLOW) {
    	for(int i = 0; i < units.size(); ++i) {
    		units[i]->sch->excess = -units[i]->perfomance_ptu / task_processing_time_expectation;
    		for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
    			units[i]->outgoing_channels[j]->scheduled_load = 0;
    		}
    	}

/*    	for(int tt = 0; tt < 30; ++tt) {
    		for(int j = 0; j < units.size(); ++j)
    		for(int i = 0; i < units[j]->sch->saturation_ratio.size(); ++i) {
    		channel* ch = units[j]->sch->saturation_ratio[i].second;
    		
			int new_f;
			if(ch->out->buffer.size() == 0) {
				new_f = ch->capacity_ptu;
			} else {
				new_f = (ch->out->sch->excess * units[j]->buffer.size()) / ch->out->buffer.size() - units[j]->sch->excess;
				if(new_f < 0) new_f = 0;
				if(new_f > ch->capacity_ptu) new_f = ch->capacity_ptu;
			}
			ch->out->sch->excess -= ch->scheduled_load / task_content_size_expectation - new_f;
			units[j]->sch->excess += ch->scheduled_load / task_content_size_expectation - new_f;
			ch->scheduled_load = new_f * task_content_size_expectation;
    		}
    	}*/
    }
    print("Milestone 2\n", 4);

    /*
	for(int i = 0; i < units.size(); ++i) {
		cout << i << ":" << units[i]->buffer.size() * task_processing_time_expectation << endl;
	}	
	*/

	if(simulation == 0) {
		return;
	}
	//print("Proceed simulation (y/n)?:", 1);
	char ch = 'y';
	//scanf("%c", &ch);
	if('A' <= ch && ch <= 'Z') ch += 'a' - 'A';
	if(ch != 'y') {
		print("Simulation aborted\n", 1);
		return;
	}
	print("Simulation started ...\n", 1);
	int current_time = 0;
	print("Simulation time: ", 4);
	print(current_time + 1, 4);
	print("\n", 4);
	print("Time output step = ", 2);
	print(time_output_step, 2);
	print("\n", 2);

	while(!tick()) {
		current_time++;
		//cerr << current_time << ": tasks resolved = " << task_resolved << endl;
		if(current_time % time_output_step == 0) {
        	print("Simulation time = ", 3);
			print(current_time, 3);
			print(": tasks resolved ", 3);
			print(task_resolved, 3);
			print("\n", 3);
		}
		print("Simulation time: ", 4);
		print(current_time + 1, 4);
		print("\n", 4);
	}
	current_time++;
	//cout << current_time << ": tasks resolved = " << task_resolved << " expected = " << task_cnt << ", time without transferring = " << m_t << endl;
	print("Simulation time = ", 2);
	print(current_time, 2);
	print(": tasks resolved ", 2);
	print(task_resolved, 2);
	print(", expected = ", 2);
	print(task_cnt, 2);
	print(", time without transferring = ", 2);
	print(m_t, 2);
	print("\n", 2);

	if(algorithm == CONSENSUS || algorithm == ADAPTIVE_FLOW) {
		print("Total scheduling time ", 2);
		print((double)scheduling_time / CLOCKS_PER_SEC, 2);
		print("s.(", 2);
		print((double)scheduling_time, 2);
		print(")\n", 2);
	} 


	int task_left_1 = 0;
    int task_left_2 = 0;
	int task_left_3 = 0;
	int task_left_4 = 0;
    for(int i = 0; i < n; ++i) {
    	task_left_1 += units[i]->buffer.size();
    	task_left_2 += units[i]->abstract_delayer->size();
    	task_left_3 += units[i]->current_task != NULL;
    	for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
    		task_left_4 += units[i]->outgoing_channels[j]->current_task != NULL;
    	}
    }

    

	//cout << "Task lost = " << task_left_1 << " + " << task_left_2 << " + " << task_left_3 << " + " << task_left_4 << endl;
	for(int i = 0; i < n; ++i) {
    	for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
    		print(i + 1, 4);
    		print("->", 4);
    		print(1 + units[i]->outgoing_channels[j]->out->id, 4);
    		print(": ", 4);
			print(units[i]->outgoing_channels[j]->transferred_load, 4);
			print(" ", 4);
			print(units[i]->outgoing_channels[j]->scheduled_load, 4);
			print("\n", 4);
    		
    		//cerr << i + 1 << "->" << 1 + units[i]->outgoing_channels[j]->out->id << ": " << 
    		//	units[i]->outgoing_channels[j]->transferred_load << " " << units[i]->outgoing_channels[j]->scheduled_load << " " << endl;
    		//cerr << "Scheduled load: " << units[i]->outgoing_channels[j]->scheduled_load << " " << graph->get_flow(n + units[i]->outgoing_channels[j]->id) << endl; 
    	}
    }


	//cout << "Current simulator time is " << current_time << endl;
}