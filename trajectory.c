#include <string.h>
#include <stdint.h>
#include <stdbool.h>
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
#include "log.h"
#include "loc.h"
#include <time.h>
#include "particle_filter.h"

void takeoff(float height);
void land(void);
void headToSetpoint(float x, float y, float z, float yaw);
void positionSet(setpoint_t *setpoint, float x, float y, float z, float yaw);
void flyCircle(point_t center, float radius, float phase, uint32_t duration_ms, uint8_t repetitions, bool clockwise);
void flyCircleStep(float radius, int step);
void flyToPoint(point_t endPos, float velocity, float yaw);
void setYaw(float yaw, float height);
void save_range(uint32_t range);

uint32_t range;
float r;
point_t p;
bool converged;

uint8_t start = 0;
uint8_t status = 0;
float PostParticles[2][NPARTICLES][2];
float sol[2] = {0, 0};
float (*p_PostParticles)[NPARTICLES][2] = PostParticles;
float (*p_sol) = sol;
int min_dist_id = 0;
float next_point[2];
float (*p_next) = next_point;

void fly_task(void* parameters) {
    DEBUG_PRINT("Fly task started!\n");
    vTaskDelay(1000);
    while(1) {

        while(!start)
            vTaskDelay(1000);
        // Reset estimator
        
        estimatorKalmanInit();
        RCC_AHB2PeriphClockCmd(RCC_AHB2Periph_RNG, ENABLE);
        RNG_Cmd(ENABLE);
        float t_start, t_end;
        t_start = (float)xTaskGetTickCount();
        
        p.z = 0.2;
        p.x = 0.0;
        p.y = 0.0;
        
        for(int k = 0; k < NPARTICLES; k++)
        {
            PostParticles[0][k][0] = MAPEDGE * ((int)RNG_GetRandomNumber()/(pow(2, 31) - 1));
            PostParticles[0][k][1] = MAPEDGE * ((int)RNG_GetRandomNumber()/(pow(2, 31) - 1));
        }
        // Take off
        status = 1;
        takeoff(0.2);
        
        status = 2;
        
        flyCircle(p, 1.5, 0.0f, 10000, 1, true);
        //flyToMinDist();
        t_end = (float)xTaskGetTickCount();
        DEBUG_PRINT("Fly to Land time: %f\n", (t_end - t_start)*portTICK_PERIOD_MS/1000);
        DEBUG_PRINT("Found Node at: x = %f, y = %f\n", sol[0], sol[1]);
        
//      Land
        status = 3;
        land();
        
        
        while(1)   vTaskDelay(1000);
    }
}

void flyCircle(point_t center, float radius, float phase, uint32_t duration_ms, uint8_t repetitions, bool clockwise){
    uint8_t step_delay_ms = 100;
    int steps = duration_ms / step_delay_ms;
    float loop_progress;
    phase = phase * (float)M_PI / 180.0f;

    // Go to start
    point_t start_point;
    start_point.x = (float)cos(phase) * radius + center.x;
    start_point.y = (float)sin(phase) * radius + center.y;
    start_point.z = center.z;

    flyToPoint(start_point, 0.3, 0.0);
    vTaskDelay(50);

    // Loop circle
    for (int j = 0; j < repetitions; j++){
        for (int i = 0; i <= steps; i++) {
            loop_progress = (float)i / (float)steps;
            float angle = 0.0f;
            if (clockwise) {
                angle = loop_progress * 2 * (float)M_PI + phase;
            } else {
                angle = - loop_progress * 2 * (float)M_PI + phase;
            }
            float x = (float)cos(angle) * radius ;//+ center.x;
            float y = (float)sin(angle) * radius ;//+ center.y;
            float z = center.z;
            
            r = range;
            
            headToSetpoint(x, y, z, 0);
            particle_filter(logGetFloat(logGetVarId("kalman", "stateX")), logGetFloat(logGetVarId("kalman", "stateY")), p_PostParticles, p_sol, r, i, steps);
            
            for(int l = 0; l < NPARTICLES; l++)
            {
                PostParticles[0][l][0] = PostParticles[1][l][0];
                PostParticles[0][l][1] = PostParticles[1][l][1];
                if(converged)
                {
                    if(MeasureDist(PostParticles[0][l][0], PostParticles[0][l][1], logGetFloat(logGetVarId("kalman", "stateX")), logGetFloat(logGetVarId("kalman", "stateY"))) - r/1000 < PostParticles[0][min_dist_id][0], PostParticles[0][min_dist_id][1], logGetFloat(logGetVarId("kalman", "stateX")), logGetFloat(logGetVarId("kalman", "stateY")) - r/1000)
                        min_dist_id = i;
                }
                
            }
            if(converged)
            {
                sol[0] = PostParticles[0][min_dist_id][0];
                sol[1] = PostParticles[0][min_dist_id][1];
                DEBUG_PRINT("Converged at step: %d\n", i);
                break;
            }
        
            vTaskDelay(step_delay_ms);
        }
    }

    p.x = sol[0];
    p.y = sol[1];
    
    flyToPoint(p, 0.7, 0);
    vTaskDelay(50);
    
}

void takeoff(float height) {
    point_t pos;
    memset(&pos, 0, sizeof(pos));
    estimatorKalmanGetEstimatedPos(&pos);

    uint32_t endheight = (uint32_t)(100 * (height - 0.3f));
    for(uint32_t i=0; i<endheight; i++) {
        headToSetpoint(pos.x, pos.y, 0.3f + (float)i / 100.0f, 0);
        estimatorKalmanGetEstimatedPos(&pos);
        vTaskDelay(30);
    }
}

void land(void) {
    point_t pos;
    memset(&pos, 0, sizeof(pos));
    estimatorKalmanGetEstimatedPos(&pos);

    float height = pos.z;
    for(int i=(int)100*height; i>18; i--) {
        headToSetpoint(pos.x, pos.y, (float)i / 100.0f, 0);
        vTaskDelay(100);
    }

    motorsSetRatio(MOTOR_M1, 0);
    motorsSetRatio(MOTOR_M2, 0);
    motorsSetRatio(MOTOR_M3, 0);
    motorsSetRatio(MOTOR_M4, 0);
    vTaskDelay(200);
}

void headToSetpoint(float x, float y, float z, float yaw) {
    setpoint_t setpoint;
    positionSet(&setpoint, x, y, z, yaw);
    commanderSetSetpoint(&setpoint, 3);
}

void positionSet(setpoint_t *setpoint, float x, float y, float z, float yaw) {
    memset(setpoint, 0, sizeof(setpoint_t));

    setpoint->mode.x = modeAbs;
    setpoint->mode.y = modeAbs;
    setpoint->mode.z = modeAbs;

    setpoint->position.x = x;
    setpoint->position.y = y;
    setpoint->position.z = z;

    setpoint->mode.yaw = modeAbs;

    setpoint->attitude.yaw = yaw;

    setpoint->mode.roll = modeDisable;
    setpoint->mode.pitch = modeDisable;
    setpoint->mode.quat = modeDisable;
}

float angle_diff(float yaw0, float yaw1) {
    float diff;
    if (yaw0 * yaw1 >= 0.0f)
        diff = yaw1 - yaw0;
    else {
        if (yaw0 >= 0.0f) {
            if (yaw0 - yaw1 > 180.0f)
                diff = 360.0f - (yaw0 - yaw1);
            else
                diff = -yaw0 + yaw1;
        }

        if (yaw0 < 0.0f) {
            if (-yaw0 + yaw1 > 180.0f)
                diff = - (360.0f - (-yaw0 + yaw1));
            else
                diff = -yaw0 + yaw1;        
        }   
    }

    return diff;
}

void setYaw(float target_yaw, float height) {
    point_t currentPos;
    memset(&currentPos, 0, sizeof(currentPos));
    estimatorKalmanGetEstimatedPos(&currentPos);

    float current_yaw = logGetFloat(logGetVarId("stateEstimate", "yaw"));
    float yaw_res = 2.0f;

    float diff = angle_diff(current_yaw, target_yaw);
    uint16_t steps = (uint16_t) ((uint16_t)(fabs(diff)) / yaw_res);

    for(uint16_t i=0; i <steps; i++) {
        current_yaw += diff / (float)steps;
        if (current_yaw > 180.0f)
            current_yaw = current_yaw - 360.0f;
        if (current_yaw < -180.0f)
            current_yaw = current_yaw + 360.0f;
        headToSetpoint(currentPos.x, currentPos.y, height, current_yaw);
        vTaskDelay(30);
    }

    for (uint16_t i = 0; i < 5; i++) {
        headToSetpoint(currentPos.x, currentPos.y, height, target_yaw);
        vTaskDelay(100);
    }
}

void flyToPoint(point_t endPos, float velocity, float yaw){
    point_t startPos;
    memset(&startPos, 0, sizeof(startPos));
    estimatorKalmanGetEstimatedPos(&startPos);

    float distX = endPos.x - startPos.x;
    float distY = endPos.y - startPos.y;
    float distance = sqrt(pow(distX, 2) + pow(distY, 2));

    float time_ms = 1000.0f * distance / velocity;
    uint16_t steps = (uint16_t) (time_ms / 50.0f);

    for (uint16_t i = 0; i < steps; i++) {
        float dx = i * distX / (float) steps;
        float dy = i * distY / (float) steps;
        headToSetpoint(startPos.x + dx, startPos.y + dy, endPos.z, yaw);
        vTaskDelay(50);
    }
    for (uint16_t i = 0; i < 5; i++) {
        headToSetpoint(endPos.x, endPos.y, endPos.z, yaw);
        vTaskDelay(100);
    }
}
void flyToMinDist()
{
    
    p.x = 0.2;
    p.y = 0.2;
    p.z = 0.2;
    
    for(int i = 0; i < NSIM; i++)
    {
        
        flyToPoint(p, 0.7, 0);
        r = range;
        
        
        particle_filter_min_dist(logGetFloat(logGetVarId("kalman", "stateX")), logGetFloat(logGetVarId("kalman", "stateY")), p_PostParticles, p_sol, r, i, p_next);
        
        if(converged)
        {
            DEBUG_PRINT("Converged at step: %d\n", i);
            p.x = next_point[0];
            p.y = next_point[1];
            break;
        }
        
        p.x = next_point[0];
        p.y = next_point[1];
        
        for(int l = 0; l < NPARTICLES; l++)
        {
            PostParticles[0][l][0] = PostParticles[1][l][0];
            PostParticles[0][l][1] = PostParticles[1][l][1];
            
        }
    }
}

//void save_range(float range){
//    r[i] = range/1000;
//    i++;
//}

PARAM_GROUP_START(CMDS)
PARAM_ADD(PARAM_UINT8, fly, &start)
PARAM_GROUP_STOP(CMDS)

LOG_GROUP_START(STATUS)
LOG_ADD(LOG_UINT8, sts, &status)
LOG_GROUP_STOP(STATUS)
