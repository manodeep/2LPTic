#pragma once

#include<stdio.h>

void initialize_powerspectrum(void);
void read_power_table(void);
double PowerSpec(double k);
void add_WDM_thermal_speeds(float *vel);
double F_Omega(double a);
double F2_Omega(double a);
double GrowthFactor(double astart, double aend);

