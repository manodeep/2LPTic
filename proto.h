#pragma once

void   print_spec(void);
int    FatalError(int errnum);
void   displacement_fields(void);
void   initialize_ffts(void);
void   set_units(void);
void   free_ffts(void);
double fnl(double x);
double periodic_wrap(double x);



