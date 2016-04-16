#include <iostream>

using namespace std;

#include "logger.h"

void print(const string msg, const int priority) {
	if(priority <= log_level) {
		cout << msg;
	}
}

void print(const char* msg, const int priority) {
	if(priority <= log_level) {
		cout << msg;
	}
}

void print(const int number, const int priority) {
	if(priority <= log_level) {
		cout << number;
	}
}

void print(const long long number, const int priority) {
	if(priority <= log_level) {
		cout << number;
	}
}


void print(const unsigned int number, const int priority) {
	if(priority <= log_level) {
		cout << number;
	}
}

void print(const Rational number, const int priority) {
	if(priority <= log_level) {
		cout << number.nom << "/" << number.denom;
	}
}

void print(const double num, const int priority) {
	if(priority <= log_level) {
		cout << num;
	}
}

void print(const char     c, const int priority) {
	if(priority <= log_level) {
		cout << c;
	}
}

void set_log_level(int level) {
	if(0 <= level && level <= 4) {
		log_level = level;
	}
}

int get_log_level() {
	return log_level;
}

