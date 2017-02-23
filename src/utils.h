#ifndef __UTILS__
#define __UTILS__

#include <time.h>









typedef void *PtrVoid;

#ifndef DEBUG
#define DEBUG 0
#endif

#define ERROR(status, ...) do {fprintf(stderr,"ERROR %s:%d in %s\n",__FILE__,__LINE__,__FUNCTION__);fprintf(stderr, __VA_ARGS__);exit(status);} while(0);
#define WARNING(...) do {time_t timer;char buffer[25];struct tm* tm_info;time(&timer);tm_info = localtime(&timer);strftime(buffer, 25, "%Y:%m:%d%H:%M:%S", tm_info);fprintf(stderr,"WARNING %s %s:%d in %s\n",buffer,__FILE__,__LINE__,__FUNCTION__);fprintf(stderr,"   ");fprintf(stderr, __VA_ARGS__);} while(0);
#define DEBUGP(...) \
  do { if (DEBUG) {time_t timer;char buffer[25];struct tm* tm_info;time(&timer);tm_info = localtime(&timer);strftime(buffer, 25, "%Y:%m:%d%H:%M:%S", tm_info);fprintf(stderr,"DEBUG %s %s:%d ",buffer, __FILE__,__LINE__);fprintf(stderr, __VA_ARGS__);} } while (0)

#endif

