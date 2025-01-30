//
//  particle_filter.h
//  
//
//  Created by Marwan Fathy on 21.10.22.
//

#ifndef particle_filter_h
#define particle_filter_h

#define SIGMA 0.2
#define K 0.3
#define K2 0.25
#define NPARTICLES 200
#define NSIM 100
#define MAPEDGE 2
#define CONV_RADIUS 0.45


void particle_filter(float Dx, float Dy, float (*PostParticles)[NPARTICLES][2], float* sol, float dist, int simstep, const int steps);
void particle_filter_min_dist(float Dx, float Dy, float (*PostParticles)[NPARTICLES][2], float* sol, float dist, int simstep, float* p_next);
float MeasureDistToDrone(float Px, float Py, float Dx, float Dy, float dist);


#endif /* particle_filter_h */
