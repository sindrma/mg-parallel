#ifndef TIMER_H__
#define TIMER_H__

#include <stdbool.h>
#define NUM_TIMERS 10

typedef struct mg_timer_s
{
    double start, end;
    bool running;
} mg_timer_t;

void init_timers();
void timer_start(int timer_id);
void timer_stop(int timer_id);
double timer_elapsed(int timer_id);

#endif
