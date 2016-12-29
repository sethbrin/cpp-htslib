#ifndef KTHREAD_H_
#define KTHREAD_H_

#ifdef __cplusplus
extern "C" {
#endif

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

#ifdef __cplusplus
}
#endif

#endif
