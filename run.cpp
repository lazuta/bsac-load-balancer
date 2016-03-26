#include <cstdlib>
#include <vector>
#include <list>
#include <ctime>
#include <queue>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include "graph.h"
#include "simulator.h"
#include "logger.h"

using namespace std;

void set_config() {
	FILE* settings = fopen("settings.cfg", "r");
	char s[100], t[50];
	while(fscanf(settings, "%s", s) == 1) {		
		int tmp;
		int equality = 0;
		while(s[equality] != 0 && s[equality] != '=') {
			equality++;
		}
		if(s[equality] == '=') {
			int i;
			for(i = 0; s[equality + i + 1] != 0; ++i) {
				t[i] = s[equality + i + 1];
			}
			s[equality] = t[i] = 0;
		} else continue;

		print(s, 4);
		print('=', 4);
		print(t, 4);
		print("\n", 4);

		if(strcmp(s, "log_level") == 0) {
			if (sscanf(t, "%d", &tmp) == 1){
				set_log_level(tmp);
			}
		} else if (strcmp(s, "task_processing_time_expectation") == 0) {
			if(sscanf(t, "%d", &tmp) == 1) {
				task_processing_time_expectation = tmp;
			}
		} else if (strcmp(s, "task_processing_time_variance") == 0) {
			if(sscanf(t, "%d", &tmp) == 1) {
				task_processing_time_variance = tmp;
			}
		} else if (strcmp(s, "task_content_size_expectation") == 0) {
			if(sscanf(t, "%d", &tmp) == 1) {
				task_content_size_expectation = tmp;
			}
		} else if (strcmp(s, "task_content_size_variance") == 0) {
			if(sscanf(t, "%d", &tmp) == 1) {
				task_content_size_variance = tmp;
			}
		}
	}	
}   

int main(int argc, char *argv[])  {
	if(argc <= 1 ) {
		return 0;
	}	
	set_config();
	print("Log level = ", 3);
	print(log_level, 3);
	print("\n", 3);

	srand(time(NULL));
	string str = "tests/";
	for(int i = 0; argv[1][i] != 0; ++i) {
		str.push_back(argv[1][i]);
	}
	simulate(str);

	return 0;
}