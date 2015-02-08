#include "includes.h"
#include "timer.h"

mg_timer_t timers[NUM_TIMERS];

inline double getTime()
{
    const double kMicro = 1.0e-6;
    struct timeval TV;

    const int RC = gettimeofday(&TV, NULL);
    if(RC == -1)
    {
        printf("ERROR: Bad call to gettimeofday\n");
        return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()


void init_timers()
{
    memset(timers, 0, sizeof(mg_timer_t)*NUM_TIMERS);
}

void timer_start(int timer_id)
{
    if (timer_id >= NUM_TIMERS) return;
    timers[timer_id].start = getTime();
    timers[timer_id].running = true;
}

void timer_stop(int timer_id)
{
    if (timer_id >= NUM_TIMERS) return;
    timers[timer_id].end = getTime();
    timers[timer_id].running = false;
}

double timer_elapsed(int timer_id)
{
    if (timer_id >= NUM_TIMERS) return -1.f;
    return (!timers[timer_id].running ? timers[timer_id].end : getTime()) - timers[timer_id].start;
}


