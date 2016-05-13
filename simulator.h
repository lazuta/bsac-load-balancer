#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <time.h>
#include <string>
#include "graph.h"
#include "utils/rational.h"


/**
 * Some forward declaration.
 */
class processing_unit;
class channel;
class delayer;
class task;
class scheduler;

/**
 * In this simulation we consider that processing time for any
 * task is a random value with uniform distribution over integer
 * values that lies in some bounded interval, that is processing
 * time will be random integer value from 
 * [expectation - variance; expectation + variance] with equal 
 * probability.
 */
static int task_processing_time_expectation = 100;
static int task_processing_time_variance = 10;
/**
 * In this simulation we consider that information that needs to be
 * transferred when giving away some task is a random value with 
 * uniform distribution over integer values that lies in some bounded 
 * interval, that is content size will be random integer value from 
 * [expectation - variance; expectation + variance] with equal 
 * probability.
 */                                           
static int task_content_size_expectation = 10;
static int task_content_size_variance = 1;

/**
 * Delayer uses this variable to determine the array size.
 */
static const int max_delay = 10;                   

/**
 * This class contains methods for emulating communication delays.
 */
class delayer {
public:	
	int cnt;
	list<task*> time_slots[max_delay];
	int current_time;

	/**
	 * Tasks to be released at current time moment.
	 */
	list<task*> get_current();

	/**
	 * Pushes delayer time one unit further and releases current tasks.
	 */
	void tick();
	
	/**
	 * Adds a task that need to be released 'delay' time units later.
	 */
	void add_task(task* t, int delay);
	
	delayer();
	
	/**
	 * Checks is delayer is empty.
	 */
    bool empty();

    /**
     * Returns the number of tasks in delayer.
     */
    int size();
};

/**
 * Class containing the main attributes of processing node.
 */
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

	processing_unit();

	/**
	 * Adds a task for this processing unit to perform.
	 */
	void add_task(task* tsk);
	
	/**
	 * Adds a communication channel to this processing unit.
	 */
    void add_channel(channel* ch);
};

/**
 * Class representing a communication channel.
 */
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

	channel();
};

/**
 * Class representing a task.
 */
class task {
public:	
	int id;
	int processing_time;
	int content_size;
	processing_unit* current_pu;
	bool resolved;

	task();
};

/**
 * Class containing basic local rules for task scheduling within a single node.
 */
class scheduler {
public:	
	vector<pair<double, channel*> > saturation_ratio;
	int task_bound;

   	bool is_saturated();
   	
   	scheduler();
};


enum balancing_algorithm_t {PARAMETRIC_FLOW, CONSENSUS};

/**
 * Simulation global variables.
 */
static int task_resolved = 0;
static double m_t = 0;
static long long task_cnt = 0;
static vector<processing_unit*> units;
static long long sum_perf = 0;
static int n = 0; // number of processing units.
static int m = 0; // number of channels.
static double tau = 0; // expected optimal time.
static balancing_algorithm_t algorithm = PARAMETRIC_FLOW;
static double alpha = 1; // consensus step size.
static clock_t scheduling_time = 0; //this variable is used to store scheduling time of algorithms (in clocks)
static int time_output_step = 1;
static int simulation = 1;
       
void read_graph(static_flow_graph* graph);

void set_consensus_step(double step);

void set_balancing_algorithm(char* algorithm_name);

void set_time_output_step(int step);

void set_simulation(int doornot);

void simulate(string path);

#endif