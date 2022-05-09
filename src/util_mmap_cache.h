#ifndef UTIL_MMAP_CACHE_H
#define UTIL_MMAP_CACHE_H

void *allocate_cache(char *variable, size_t length);
void deallocate_cache(void *addr, size_t length);

#endif  // UTIL_MMAP_CACHE_H
