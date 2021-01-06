/*
   based on: https://raw.githubusercontent.com/aramadia/udacity-cs344/master/Unit2%20Code%20Snippets/gputimer.h
*/

#ifndef __GPU_TIMER_H__
#define __GPU_TIMER_H__

#include <cuda.h>

class GpuTimer
{
public:
      GpuTimer()
      {
            cudaEventCreate(&start_);
            cudaEventCreate(&stop_);
      }
 
      ~GpuTimer()
      {
            cudaEventDestroy(start_);
            cudaEventDestroy(stop_);
      }
 
      void start()
      {
            cudaEventRecord(start_, 0);
      }
 
      void stop()
      {
            cudaEventRecord(stop_, 0);
      }
 
      float elapsed()
      {
            float elapsed;
            cudaEventSynchronize(stop_);
            cudaEventElapsedTime(&elapsed, start_, stop_);
            return elapsed;
      }

private:
      cudaEvent_t start_;
      cudaEvent_t stop_;
};

#endif  /* __GPU_TIMER_H__ */
