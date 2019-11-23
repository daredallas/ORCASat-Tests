/*
* Cubesat_Mechanics.c
* Cubesat mechanics calculations
*/

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "./Reference_Libraries/Matrix_Math_Library/matrix.h"
#include "./Reference_Libraries/Matrix_Math_Library/Quaternion.h"


typedef struct Cube_Mechs{
  gsl_vector* acceleration;
  gsl_vector* velocity;
  gsl_vector* position;
  Quaternion q_BI;
  gsl_vector* V_b;
  gsl_matrix* DCM_BI;
  gsl_vector* omega_b;
  gsl_vector* angular_momentum;
  gsl_vector* next_integrator;
  Quaternion prev_q_BI;
}Cube_Mechs;


/*
* Applies Runge-Kutta 4th order formula
*/
double rk4(double timestamp, double value){
  return (5*timestamp*timestamp - value) / exp(timestamp+value);
}


/*
* Executes the runge-kutta 4th order ODE on the given vector with number of elements provided by the elements parameter
* Inputs: Previous acceleration vector, timestamp, and number of elements in acceleration vector
* Outputs: A new acceleration vector
*/
gsl_vector* runge_kutta4(gsl_vector* acceleration, double timestamp, int elements){

  double step_size = .1;
  double F1, F2, F3, F4;
  double values[elements];

  for(int i=0; i<elements; i++){
    double cur = gsl_vector_get(acceleration,i);

    F1 = .1*rk4(timestamp, cur);
    F2 = .1*rk4((timestamp + step_size/2), cur + F1/2);
    F3 = .1*rk4((timestamp + step_size/2), cur + F2/2);
    F4 = .1*rk4((timestamp + step_size), cur + F3);

    values[i] = cur + ((F1+F2+F3+F4)/6);
  }

  return init_vector_custom(elements, values);
}


/*
* Performs all force-based calculations for the cubesat mechanics block in the Simulink model
* Inputs: Current force, timestamp, and struct used to wrap all data
* Outputs: Acceleration, velocity, and position, wrapped in the Cube_Mechs struct
*/
void force_calculations(gsl_vector* force, double timestamp, Cube_Mechs *cm){
  double mass = 3.6;
  gsl_vector* normalized_force = scalar_multiplication_vector(force, (1/mass));
  cm->acceleration = normalized_force;
  cm->velocity = runge_kutta4(normalized_force, timestamp, 3);
  cm->position = runge_kutta4(cm->velocity, timestamp, 3);
}


/*
* Executes the creation of the subsystem2 vector from the cubesat mechanics block in the Simulink model
* Inputs: Integrator vector and wheel moment vecotr
* Outputs: subsystem2 vector
*/
gsl_vector* subsystem2(gsl_vector* integrator, gsl_vector* wheel_moment_new){
  double values[3] = {(gsl_vector_get(integrator, 1)*gsl_vector_get(wheel_moment_new, 2))-(gsl_vector_get(integrator, 2)*gsl_vector_get(wheel_moment_new, 1)),
                      (gsl_vector_get(integrator, 2)*gsl_vector_get(wheel_moment_new, 0))-(gsl_vector_get(integrator, 0)*gsl_vector_get(wheel_moment_new, 2)),
                      (gsl_vector_get(integrator, 0)*gsl_vector_get(wheel_moment_new, 1))-(gsl_vector_get(integrator, 1)*gsl_vector_get(wheel_moment_new, 0))};
  return init_vector_custom(3, values);
}


/*
* Performs all force-based calculations for the cubesat mechanics block in the Simulink model
* Inputs: Current force, timestamp, and struct used to wrap all data
* Outputs: Acceleration, velocity, and position, wrapped in the Cube_Mechs struct
*/
Quaternion quaternion_kinematics(gsl_vector* omega, Quaternion prev_q_BI, double timestamp){
  gsl_matrix* skew = skew_symmetric_vector(omega);
  skew = scalar_multiplication_matrix(skew, -1);

  gsl_matrix* omega_vertical = vector_to_matrix(3, omega);
  omega_vertical = concatonate_matrices_horizontal(skew, omega_vertical);

  double omega_horizontal_values[1][4] = {gsl_vector_get(omega,0), gsl_vector_get(omega,1), gsl_vector_get(omega,2), 0};
  gsl_matrix* omega_horizontal = init_matrix_custom(1, 4, omega_horizontal_values);

  gsl_matrix* Omega = concatonate_matrices_vertical(omega_vertical, omega_horizontal);

  gsl_vector* MM = product_vector_matrix(Omega, prev_q_BI.values);
  MM = scalar_multiplication_vector(MM, .5);

  Quaternion q;
  q.values = normalize_vector(runge_kutta4(MM, timestamp, 4));
  return q;
}


/*
* Performs all rotation-based calculations for the cubesat mechanics block in the Simulink model
* Inputs: Current moment, control Torque, momentum wheel torque, momentum wheel momentum, the previous integration value from the runge-kutta ODE, the previous Body->ECI quaternion, timestamp, and Cube_Mechs struct
* Outputs: Body->ECI quaternion, Body->ECI DCM, velocity in the body frame, omega in the body frame, angular momentum, and the next integration value for use in the runge-kutta calculation
*/
void rotation_calculations(gsl_vector* moment, gsl_vector* control_torque, gsl_vector* wheel_torque, gsl_vector* wheel_moment, gsl_vector* integrator, Quaternion prev_q_BI, double timestamp, Cube_Mechs *cm){
  integrator = runge_kutta4(integrator, timestamp, 4);

  double j_values[3][3] = {{.008, 0, 0}, {0, .008, 0}, {0, 0, .003}};
  gsl_matrix* J_0 = init_matrix_custom(3, 3, j_values);

  gsl_vector* MM1 = product_vector_matrix(J_0, integrator);

  gsl_vector* wheel_moment_new = add_vector(wheel_moment, MM1);

  gsl_vector* C = subsystem2(integrator, wheel_moment_new);

  gsl_vector* dynamics = add_vector(moment, control_torque);
  dynamics = subtract_vector(dynamics, wheel_torque);
  dynamics = subtract_vector(wheel_moment, C);

  gsl_vector* new_integrator = product_vector_matrix(invert_matrix(J_0), dynamics);

  cm->q_BI = quaternion_kinematics(integrator, prev_q_BI, timestamp);
  cm->DCM_BI = quaternion_to_dcm(cm->q_BI);
  cm->V_b = product_vector_matrix(cm->DCM_BI, cm->velocity);

  cm->omega_b = integrator; //omega
  cm->angular_momentum = C; //Angular momentum
  cm->next_integrator = new_integrator; //Next 1/s value
}


void Cubesat_Mechanics(gsl_vector* force, gsl_vector* moment, gsl_vector* control_torque, gsl_vector* wheel_torque, gsl_vector* wheel_moment, gsl_vector* integrator, Quaternion prev_q_BI, double timestamp){
  Cube_Mechs cm;
  force_calculations(force, timestamp, &cm);
  rotation_calculations(moment, control_torque, wheel_torque, wheel_moment, integrator, prev_q_BI, timestamp, &cm);
}
