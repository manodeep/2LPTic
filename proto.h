#pragma once


#define ASSERT_ALLOC(cond)												\
  do {																	\
	if (!(cond)) {														\
	  printf("failed to allocate %g Mbyte on Task %d \n", bytes / (1024.0 * 1024.0), ThisTask); \
	  printf("Error in file: %s\tfunc: %s\tline: %d with expression `"#cond"'\n", __FILE__, __FUNCTION__, __LINE__); \
	  printf("bailing out.\n");											\
	  FatalError(1);													\
	} else {															\
	  if(ThisTask == 0)													\
		printf("\nallocated %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask); \
	}																	\
  } while (0)															\


#if 0
#define ASSERT_ALLOC(cond) {                                                                                  \
   if(cond)                                                                                                   \
    {                                                                                                         \
      if(ThisTask == 0)                                                                                       \
	printf("\nallocated %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                     \
    }                                                                                                         \
  else                                                                                                        \
    {                                                                                                         \
      printf("failed to allocate %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                \
      printf("bailing out.\n");                                                                               \
      FatalError(1);                                                                                          \
    }                                                                                                         \
}
#endif

void   print_spec(void);
int    FatalError(int errnum);
void   displacement_fields(void);
void   initialize_ffts(void);
void   set_units(void);
void   free_ffts(void);
double fnl(double x);
double periodic_wrap(double x);



