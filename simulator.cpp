#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <list>
#include <conio.h>

using namespace std;

#include "simulator.h"
#include "logger.h"


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

void read_graph(static_flow_graph* graph) {
	int s, t;
	scanf("%d%d", &n, &m);
	s = n, t = n + 1;
	task_cnt = 0;
	sum_perf = 0;
	
	m_t = 0;

	for(int i = 0; i < n; ++i) {
		int p, q;
		scanf("%d%d", &p, &q);
		processing_unit* u = new processing_unit();
		u->sch = new scheduler();
		u->id = i;
		u->abstract_delayer = new delayer();
		units.push_back(u);
		u->perfomance_ptu = p * task_processing_time_expectation;
		sum_perf += p;
		double load = 0;
		for(int j = 0; j < q; ++j) {
			task* tsk = new task();
			tsk->current_pu = u;
			tsk->id = task_cnt++;
			u->add_task(tsk);
			load += tsk->processing_time;
		}
		m_t = max(m_t, load / u->perfomance_ptu);
	}
	for(int i = 0; i < m; ++i) {
		int q, w; double c;
		scanf("%d%d%lf", &q, &w, &c);
		
		channel* ch = new channel();
		ch->in = units[q - 1];
		ch->out = units[w - 1];
		ch->id = i;
		ch->delay = 1;
		ch->capacity_ptu = c * task_content_size_expectation;

		units[q - 1]->add_channel(ch);
	}	

	graph->initialize_graph(n + 2, 2 * n + m, s, t, 1e-6);
	
	for(int i = 0; i < n; ++i) {
		graph->add_edge(i, t, units[i]->perfomance_ptu / task_processing_time_expectation, 0);
    }

    for(int i = 0; i < n; ++i) {
    	for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
		    graph->add_edge(i, units[i]->outgoing_channels[j]->out->id, units[i]->outgoing_channels[j]->capacity_ptu / task_content_size_expectation, 0);
    	}
    }

    for(int i = 0; i < n; ++i) {
		graph->add_edge(s, i, 0, units[i]->buffer.size());
    }

    //cout << "Flow calculation" << endl;
    print("Flow calculation\n", 2);
    clock_t cur_time = clock();
    tau = 1.0 / graph->leftmost_breakpoint((double)sum_perf / task_cnt);
    //cout << "Done" << endl;
    print("Done\n", 2);
    //cout << "Calculation of a control took " << (double)(clock() - cur_time) / CLOCKS_PER_SEC << "s." << endl;
	print("Flow calculation took ", 2);
	print((double)(clock() - cur_time) / CLOCKS_PER_SEC, 2);
	print("Number of pushes: ", 2);
	print(graph->pushes(), 2);
	print("\n", 2);
	print("Nubmer of relabels: ", 2);
	print(graph->relabels(), 2);
	print("\n", 2);
	print("s.\n", 2);
	//cerr << "Calculation of a control took " << (double)(clock() - cur_time) / CLOCKS_PER_SEC << "s." << endl;
	
	//cout << "Found tau = " << tau << endl;
	print("Found tau = ", 2);
	print(tau, 2);
	print("\n", 2);

	graph->show(4);
}

void tick(processing_unit* unit) {
	unit->current_progress += unit->perfomance_ptu;
	if(unit->current_task != NULL && unit->current_progress >= unit->current_task->processing_time) {
		unit->current_progress -= unit->current_task->processing_time;
		unit->current_task->resolved = true;
		task_resolved++;
		//cerr << "Task #" << unit->current_task->id << " is resolved!!! Congratulations!!!" << endl;
		unit->current_task = NULL;
	}
	if(unit->current_task == NULL) {
		while(!unit->buffer.empty() && unit->current_progress >= unit->buffer.front()->processing_time) {
			unit->current_progress -= unit->buffer.front()->processing_time;
			unit->buffer.front()->resolved = true;
			task_resolved++;
			//cerr << "Task #" << unit->buffer.front()->id << " is resolved!!! Congratulations!!!" << endl;
			unit->buffer.pop();    				
		}
		if(unit->buffer.empty()) {
			unit->current_progress = 0;
		} else {
			unit->current_task = unit->buffer.front();
			unit->buffer.pop();
		}
	}                   	
	
	
	for(int i = 0; i < unit->sch->saturation_ratio.size(); ++i) {
		channel* ch = unit->sch->saturation_ratio[i].second;
		ch->current_progress += ch->capacity_ptu;
		if(ch->current_task != NULL && ch->current_progress >= ch->current_task->content_size) {
			ch->current_progress -= ch->current_task->content_size;
			ch->out->abstract_delayer->add_task(ch->current_task, ch->delay);
			ch->current_task = NULL;
		}                                       
		if(ch->current_task == NULL) {
			while(!unit->buffer.empty() && ch->transferred_load < ch->scheduled_load && ch->current_progress >= unit->buffer.front()->content_size) {
				ch->current_progress -= unit->buffer.front()->content_size;
				ch->transferred_load += unit->buffer.front()->content_size;
				ch->out->abstract_delayer->add_task(unit->buffer.front(), ch->delay);
				unit->buffer.pop();
			}
			if(unit->buffer.empty() || ch->transferred_load >= ch->scheduled_load) {
				ch->current_progress = 0;
			} else {
				ch->current_task = unit->buffer.front();
				ch->transferred_load += unit->buffer.front()->content_size;
				unit->buffer.pop();
			}
			if(ch->scheduled_load == 0) {
				unit->sch->saturation_ratio[i].first = 1;
			} else {
				unit->sch->saturation_ratio[i].first = (double)ch->transferred_load / ch->scheduled_load;        		
			}
		}
	}
	sort(unit->sch->saturation_ratio.begin(), unit->sch->saturation_ratio.end());	
	list<task*> tmp = unit->abstract_delayer->get_current();

	while(!tmp.empty()) {
		unit->buffer.push(tmp.front());
		//cerr << "Task #" << tmp.front()->id << " is transferred to unit #" << unit->id << endl;
		tmp.pop_front();
	}
	unit->abstract_delayer->tick();

	unit->status = ((int)(!unit->buffer.empty() || (unit->current_task != NULL) || (!unit->abstract_delayer->empty())) << 1) + (!unit->sch->is_saturated());
	//cerr << "Unit #" << unit->id << " status is " << unit->status << endl;
}


bool tick() {
	bool done = true;
	for(int i = 0; i < units.size(); ++i) {
		tick(units[i]);
		if(units[i]->status > 1) done = false;
	}
	return done;	
}


void simulate(string path) {
	print(path, 4);
	print("\n", 4);	
	if( ! (freopen(path.c_str(), "r", stdin)) ) {
		//cout << "Test doesn't exist" << endl;
		print("Test doesn't exist\n", 0);		
		return;
	}	

	static_flow_graph* graph;
	graph = new static_flow_graph();
    read_graph(graph);

	print("Milestone 1\n", 4);

	print("Number of processing units: ", 4);
	print(units.size(), 4);
	print("\n", 4);

	print("Number of communication channels: ", 4);
	print(m, 4);
	print("\n", 4);

    //graph->show(4);

	
	int num = 0;
    for(int i = 0; i < n; ++i) {
    	for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
    		print("Trying to set flow on (", 4);
    		print(units[i]->id, 4);
    		print(", ", 4);
    		print(units[i]->outgoing_channels[j]->id, 4);
			print(") = ", 4);
			print(floor(graph->get_flow(n + num) * tau + 0.5) * task_content_size_expectation, 4);
			print("\n", 4);
			print("Success!!!\n", 4);
    		units[i]->outgoing_channels[j]->scheduled_load = floor(graph->get_flow(n + num) * tau + 0.5) * task_content_size_expectation;
    		//cerr << "Scheduled load: " << units[i]->outgoing_channels[j]->scheduled_load << " " << graph->get_flow(n + num) * tau << endl;
    		num++; 
    	}
    }       

    print("Milestone 2\n", 4);

    /*
	for(int i = 0; i < units.size(); ++i) {
		cout << i << ":" << units[i]->buffer.size() * task_processing_time_expectation << endl;
	}	
	*/

	print("Press any key to begin simulation\n", 2);
	_getch();
	int current_time = 0;
	while(!tick()) {
		current_time++;
		//cerr << current_time << ": tasks resolved = " << task_resolved << endl;
		print("Simulation time = ", 3);
		print(current_time, 3);
		print(": tasks resolved ", 3);
		print(task_resolved, 3);
		print("\n", 3);
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