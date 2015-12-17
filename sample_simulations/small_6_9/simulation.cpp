#include <cstdlib>
#include <vector>
#include <list>
#include <ctime>
#include <queue>
#include <iostream>
#include <algorithm>
#include <cstdio>

using namespace std;

#include "flow_graph.cpp"

/**
 * ptu stands for "per time unit".
 */

/**
 * In this simulation we consider that processing time for any
 * task is a random value with uniform distribution over integer
 * values that lies in some bounded interval, that is processing
 * time will be random integer value from 
 * [expectation - variance; expectation + variance] with equal 
 * probability.
 */
const int task_processing_time_expectation = 100;
const int task_processing_time_variance = 10;
/**
 * In this simulation we consider that information that needs to be
 * transferred when giving away some task is a random value with 
 * uniform distribution over integer values that lies in some bounded 
 * interval, that is content size will be random integer value from 
 * [expectation - variance; expectation + variance] with equal 
 * probability.
 */                                           
const int task_content_size_expectation = 10;
const int task_content_size_variance = 1;

/**
 * Maximum delay for communication channels.
 * For this simulation we consider random integer values
 * of delay in range [1; max_delay].
 */
const int max_delay = 10;

class processing_unit;
class channel;
class delayer;
class task;
class scheduler;

class channel {
public:
	int id;
	processing_unit* in;
	processing_unit* out;
	int current_progress;
	int capacity_ptu;
	int delay;
	int scheduled_load;
	int transferred_load;
	task* current_task;

	channel() {
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
};

class task {
public:	
	int id;
	int processing_time;
	int content_size;
	processing_unit* current_pu;
	bool resolved;

	task() {
		id = 0;
		processing_time = task_processing_time_expectation + task_processing_time_variance - (rand() % (2 * task_processing_time_variance + 1));
		content_size = task_content_size_expectation + task_content_size_variance - (rand() % (2 * task_content_size_variance + 1));
		current_pu = NULL;
		resolved = false;
	}
};

class delayer {
public:	
	int cnt;
	list<task*> time_slots[max_delay];
	int current_time;
	list<task*> get_current() {
		return time_slots[current_time];
	}
	void tick() {
		cnt -= time_slots[current_time].size();
		time_slots[current_time].clear();
		current_time = (current_time + 1) % max_delay;
	}
	void add_task(task* t, int delay) {
		cnt++;
		time_slots[(current_time + delay) % max_delay].push_back(t);
	}

    delayer() {
    	current_time = 0;
    	cnt = 0;
    }

    bool empty() {
    	return cnt == 0;    	
    }                               
};

class scheduler {
public:	
	vector<pair<double, channel*> > saturation_ratio;
	int task_bound;

   	bool is_saturated() {
   		for(int i = 0; i < saturation_ratio.size(); ++i) {
   			if(saturation_ratio[i].second->scheduled_load < saturation_ratio[i].second->transferred_load ||
   				saturation_ratio[i].second->current_task != NULL) {
   				return false;
   			}
   		}
   		return true;
   	}

   	scheduler() {
   		task_bound = -1;
   	}
};

class processing_unit {
public:	
	int id;
	int perfomance_ptu;
	vector<channel*> outgoing_channels;
	queue<task*> buffer;
	task* current_task;
	scheduler* sch;
	delayer* abstract_delayer;

	int current_progress;
	/**
	 * 0 - buffer is empty, channels are saturated.
	 * 1 - buffer is empty, channels are not saturated.
	 * 2 - buffer is not empty, channels are saturated.
	 * 3 - buffer is not empty, channels are not saturated.
	 *
	 * Basicly, job package is done if every processing unit has 0 or 1 state.
	 * Ideally, every processing unit should have 0 state at the end of the work,
	 * but this is not critical as channels might be not completely saturated due
	 * to uncertainties of the model.  
	 */
	int status;

	processing_unit() {
		id = 0;
		perfomance_ptu = 0;
		current_task = NULL;
		sch = NULL;
		abstract_delayer = NULL;
	}

	void add_task(task* tsk) {
		buffer.push(tsk);
	}

	void add_channel(channel* ch) {
		outgoing_channels.push_back(ch);
		sch->saturation_ratio.push_back(make_pair(1, ch));
	}
};

int task_resolved = 0;

void tick(processing_unit* unit) {
	unit->current_progress += unit->perfomance_ptu;
	if(unit->current_task != NULL && unit->current_progress >= unit->current_task->processing_time) {
		unit->current_progress -= unit->current_task->processing_time;
		unit->current_task->resolved = true;
		task_resolved++;
		cerr << "Task #" << unit->current_task->id << " is resolved!!! Congratulations!!!" << endl;
		unit->current_task = NULL;
	}
	if(unit->current_task == NULL) {
		while(!unit->buffer.empty() && unit->current_progress >= unit->buffer.front()->processing_time) {
			unit->current_progress -= unit->buffer.front()->processing_time;
			unit->buffer.front()->resolved = true;
			task_resolved++;
			cerr << "Task #" << unit->buffer.front()->id << " is resolved!!! Congratulations!!!" << endl;
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
		cerr << "Task #" << tmp.front()->id << " is transferred to unit #" << unit->id << endl;
		tmp.pop_front();
	}
	unit->abstract_delayer->tick();

	unit->status = ((int)(!unit->buffer.empty() || (unit->current_task != NULL) || (!unit->abstract_delayer->empty())) << 1) + (!unit->sch->is_saturated());
	cerr << "Unit #" << unit->id << " status is " << unit->status << endl;
}

vector<processing_unit*> units;

bool tick() {
	bool done = true;
	for(int i = 0; i < units.size(); ++i) {
		tick(units[i]);
		if(units[i]->status > 1) done = false;
	}
	return done;	
}

int main() {
	freopen(".in", "rt", stdin);
	freopen(".out", "wt", stdout);

	srand(time(NULL));

	int n, m, s, t;
	scanf("%d%d", &n, &m);
	s = n, t = n + 1;
	int task_cnt = 0;
	int sum_perf = 0;
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
		for(int j = 0; j < q; ++j) {
			task* tsk = new task();
			tsk->current_pu = u;
			tsk->id = task_cnt++;
			u->add_task(tsk);
			cerr << "Task #" << task_cnt << ":" << tsk->content_size << " " << tsk->processing_time << endl;
		}
		//add_edge(i, t, perf[i]);
		//add_edge(t, i, 0);
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
		
		//add_edge(q - 1, w - 1, c);
		//add_edge(w - 1, q - 1, 0);
	}
	
	for(int i = 0; i < n; ++i) {
		//add_edge(s, i, query[i]);
		//add_edge(i, s, 0);		
	}

    static_flow_graph* graph = new static_flow_graph();

    graph->initialize_graph(n + 2, 2 * n + m, s, t, 1e-6);
    for(int i = 0; i < n; ++i) {
    	graph->add_edge(i, t, 0, units[i]->perfomance_ptu / task_processing_time_expectation);
    }

    for(int i = 0; i < n; ++i) {
    	for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
    		graph->add_edge(i, units[i]->outgoing_channels[j]->out->id, 0, units[i]->outgoing_channels[j]->capacity_ptu / task_content_size_expectation);
    	}
    }

    for(int i = 0; i < n; ++i) {
    	graph->add_edge(s, i, units[i]->buffer.size(), 0);
    }

    double tau = graph->rightmost_breakpoint((double)task_cnt / sum_perf);
	
	cerr << "Found tau = " << tau << endl;

    int num = 0;
    for(int i = 0; i < n; ++i) {
    	for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
    		units[i]->outgoing_channels[j]->scheduled_load = floor(graph->get_flow(n + num++) + 0.5) * task_content_size_expectation;
    		cerr << "Scheduled load: " << units[i]->outgoing_channels[j]->scheduled_load << " " << graph->get_flow(n + units[i]->outgoing_channels[j]->id) << endl; 
    	}
    }


	for(int i = 0; i < units.size(); ++i) {
		cout << i << ":" << units[i]->buffer.size() * task_processing_time_expectation << endl;
	}	

	int current_time = 0;
	while(!tick()) {
		current_time++;
		cerr << current_time << ": tasks resolved = " << task_resolved << endl;
	}
	cerr << current_time << ": tasks resolved = " << task_resolved << endl;
	
	for(int i = 0; i < n; ++i) {
    	for(int j = 0; j < units[i]->outgoing_channels.size(); ++j) {
    		cerr << i + 1 << "->" << 1 + units[i]->outgoing_channels[j]->out->id << ": " << units[i]->outgoing_channels[j]->transferred_load << " " << units[i]->outgoing_channels[j]->scheduled_load << " " << endl;
    		//cerr << "Scheduled load: " << units[i]->outgoing_channels[j]->scheduled_load << " " << graph->get_flow(n + units[i]->outgoing_channels[j]->id) << endl; 
    	}
    }


	cout << current_time << endl;
}