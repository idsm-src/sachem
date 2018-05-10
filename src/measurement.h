#ifndef MEASUREMENT_H_
#define MEASUREMENT_H_

#include <stdint.h>
#include <sys/time.h>


static inline struct timeval time_get()
{
    struct timeval value;
    gettimeofday(&value, NULL);
    return value;
}


static inline int64_t time_spent(struct timeval begin, struct timeval end)
{
    return ((int64_t) end.tv_sec - (int64_t) begin.tv_sec) * 1000000 + ((int64_t) end.tv_usec - (int64_t) begin.tv_usec);
}


static inline double time_to_ms(int64_t time)
{
    return time / 1000.0;
}

#endif /* MEASUREMENT_H_ */
