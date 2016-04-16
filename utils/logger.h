#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include "rational.h"

/**
 * 0  :   No log.
 * 1-3: Some log.
 * 4  : full log.
 * Use priority parameter of print with respect to these values!
 */
static int log_level = 2;

void print(const string msg, const int priority);
void print(const char*  msg, const int priority);
void print(const int number, const int priority);
void print(const long long number, const int priority);
void print(const unsigned int number, const int priority);
void print(const Rational number, const int priority);
void print(const double num, const int priority);
void print(const char     c, const int priority);     

void set_log_level(int level);
int get_log_level();

#endif