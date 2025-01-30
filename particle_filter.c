//
//  particle_filter.c
//  
//
//  Created by Marwan Fathy on 21.10.22.
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "particle_filter.h"
#include "debug.h"
#include "log.h"
#include <time.h>
#include "deck.h"
#include "debug.h"
#include <string.h>
#include <math.h>
#include "system.h"
#include "FreeRTOS.h"
#include "task.h"
#include "motors.h"
#include "estimator_kalman.h"
#include "commander.h"
#include "param.h"
#include "uwb_api.h"
#include "loc.h"
#include "../../../crazyflie-firmware/src/lib/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_rng.c"

extern bool converged = false;

float MeasureDist(float Px, float Py, float Dx, float Dy)
{
    return sqrt(pow(Px - Dx, 2) + pow(Py - Dy, 2));
}

void particle_filter(float Dx, float Dy, float (*PostParticles)[NPARTICLES][2], float* sol, float dist, int simstep, const int steps)
{
    float beta = 0, betasCumSum[NPARTICLES] = {0};
    float sum = 0;
    float random = 0;
    float max_x = 0, min_x = 0, max_y = 0, min_y = 0, avg_x = 0, avg_y = 0;
    int temp = 0;
    float StdRough[2] = {0};
    
    for( int i = 0; i < NPARTICLES; i++)
    {
        beta = exp(- pow(dist/1000 - MeasureDist(*(*(*(PostParticles) + i)), *(*(*(PostParticles) + i) + 1), Dx, Dy), 2)/(2 * pow(SIGMA, 2)));
        if(MeasureDist(*(*(*(PostParticles) + i)), *(*(*(PostParticles) + i) + 1), Dx, Dy) - dist/1000 <
           MeasureDist(*(*(*(PostParticles) + temp)), *(*(*(PostParticles) + temp) + 1), Dx, Dy) - dist/1000)
            temp = i;
        sum += beta;
        betasCumSum[i] = sum;
    }
    
    for(int i = 0; i < NPARTICLES; i++)
        betasCumSum[i] = betasCumSum[i]/sum;
    
    for(int j = 0; j < NPARTICLES; j++)
    {
        random = 0.5 * ((int)RNG_GetRandomNumber()/(pow(2, 31) - 1)) + 0.5;
        
        for(int id = 0; id < NPARTICLES; id ++)
            if(betasCumSum[id] >= random)
            {
                temp = id;
                break;
            }
        
        *(*(*(PostParticles + 1) + j)) = *(*(*(PostParticles) + temp));
        *(*(*(PostParticles + 1) + j) + 1)= *(*(*(PostParticles) + temp) + 1);
        
        if(j == 0)
        {
            max_x = *(*(*(PostParticles + 1)));
            max_y = *(*(*(PostParticles + 1)) + 1);
            min_x = max_x;
            min_y = max_y;
        }
        else
        {
            if(*(*(*(PostParticles + 1) + j)) > max_x)
                max_x = *(*(*(PostParticles + 1) + j));
            
            if(*(*(*(PostParticles + 1) + j)) < min_x)
                min_x = *(*(*(PostParticles + 1) + j));
            
            if(*(*(*(PostParticles + 1) + j) + 1) > max_y)
                max_y = *(*(*(PostParticles + 1) + j) + 1);
            
            if(*(*(*(PostParticles + 1) + j) + 1) < min_y)
                min_y = *(*(*(PostParticles + 1) + j) + 1);
        }
        
    }
    
    if(max_x - min_x < CONV_RADIUS && max_y - min_y < CONV_RADIUS)
        converged = true;
    
    StdRough[0] = K * (max_x - min_x) * 1/sqrt(NPARTICLES);
    StdRough[1] = K * (max_y - min_y) * 1/sqrt(NPARTICLES);
    
    
    for(int j = 0; j < NPARTICLES; j++)
    {
        
        *(*(*(PostParticles + 1) + j)) += StdRough[0] * 3 * ((int)RNG_GetRandomNumber()/(pow(2, 31) - 1));
        *(*(*(PostParticles + 1) + j) + 1) += StdRough[1] * 3 * ((int)RNG_GetRandomNumber()/(pow(2, 31) - 1));
    }
    
    if(simstep == steps)
    {
        for(int j = 0; j < NPARTICLES; j ++)
        {
            avg_x += *(*(*(PostParticles + 1) + j));
            avg_y += *(*(*(PostParticles + 1) + j) + 1);
        }
        
        *sol = avg_x;
        *(sol + 1) = avg_y;
    }
    
}

void particle_filter_min_dist(float Dx, float Dy, float (*PostParticles)[NPARTICLES][2], float* sol, float dist, int simstep, float* p_next)
{
    float beta = 0, betasCumSum[NPARTICLES] = {0};
    float sum = 0;
    float random = 0;
    float max_x = 0, min_x = 0, max_y = 0, min_y = 0, avg_x = 0, avg_y = 0;
    int temp = 0;
    float StdRough[2] = {0};
    
    for( int i = 0; i < NPARTICLES; i++)
    {
        beta = exp(- pow(dist/1000 - MeasureDist(*(*(*(PostParticles) + i)), *(*(*(PostParticles) + i) + 1), Dx, Dy), 2)/(2 * pow(SIGMA, 2)));
        if(MeasureDist(*(*(*(PostParticles) + i)), *(*(*(PostParticles) + i) + 1), Dx, Dy) - dist/1000 <
           MeasureDist(*(*(*(PostParticles) + temp)), *(*(*(PostParticles) + temp) + 1), Dx, Dy) - dist/1000)
            temp = i;
        sum += beta;
        betasCumSum[i] = sum;
    }
    
    *p_next = *(*(*(PostParticles) + temp));
    *(p_next + 1) = *(*(*(PostParticles) + temp) + 1);
    temp = 0;
    
    for(int i = 0; i < NPARTICLES; i++)
        betasCumSum[i] = betasCumSum[i]/sum;
    
    
    for(int j = 0; j < NPARTICLES; j++)
    {
        random = 0.5 * ((int)RNG_GetRandomNumber()/(pow(2, 31) - 1)) + 0.5;
        
        for(int id = 0; id < NPARTICLES; id ++)
            if(betasCumSum[id] >= random)
            {
                temp = id;
                break;
            }
        
        *(*(*(PostParticles + 1) + j)) = *(*(*(PostParticles) + temp));
        *(*(*(PostParticles + 1) + j) + 1)= *(*(*(PostParticles) + temp) + 1);
        
        if(j == 0)
        {
            max_x = *(*(*(PostParticles + 1)));
            max_y = *(*(*(PostParticles + 1)) + 1);
            min_x = max_x;
            min_y = max_y;
        }
        else
        {
            if(*(*(*(PostParticles + 1) + j)) > max_x)
                max_x = *(*(*(PostParticles + 1) + j));
            
            if(*(*(*(PostParticles + 1) + j)) < min_x)
                min_x = *(*(*(PostParticles + 1) + j));
            
            if(*(*(*(PostParticles + 1) + j) + 1) > max_y)
                max_y = *(*(*(PostParticles + 1) + j) + 1);
            
            if(*(*(*(PostParticles + 1) + j) + 1) < min_y)
                min_y = *(*(*(PostParticles + 1) + j) + 1);
        }
        
    }
    
    if(max_x - min_x < CONV_RADIUS && max_y - min_y < CONV_RADIUS)
        converged = true;
    
    StdRough[0] = K2 * (max_x - min_x) * 1/sqrt(NPARTICLES);
    StdRough[1] = K2 * (max_y - min_y) * 1/sqrt(NPARTICLES);
    
    
    for(int j = 0; j < NPARTICLES; j++)
    {
        
        *(*(*(PostParticles + 1) + j)) += StdRough[0] * 3 * ((int)RNG_GetRandomNumber()/(pow(2, 31) - 1));
        *(*(*(PostParticles + 1) + j) + 1) += StdRough[1] * 3 * ((int)RNG_GetRandomNumber()/(pow(2, 31) - 1));
    }
    *sol = *p_next;
    *(sol + 1) = *(p_next + 1);
    
}
