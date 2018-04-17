#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

//------------------------------------------------------------------//
//-----------------------------FUNCTIONS----------------------------//
//------------------------------------------------------------------//


double gaussian(double x);
//double Gaussian(double x);
//double MaxwBoltz(double x);


//------------------------------------------------------------------//
//-------------------------------MAIN-------------------------------//
//------------------------------------------------------------------//


int main ()
{

//------------------------------------------------------------------//
//----------------------------GLOBALS-------------------------------//
//------------------------------------------------------------------//

	double Secondary_I_Mass = 1.674e-27; // secondary ion mass [kilograms]
	double e_Charge=1.602e-19; // charge state [Coulombs]
	double t_Init_Pulse = 97e4;
	double Extraction_Voltage=2000; // extraction voltage [Volts]
	double Colimator_d = 5e-3; // maximum distance [metres]
	double t_Step=1e-10; //time step[s]
	double Central_Velocity = 100; //Medium Initial Velocity [m/s]
	double Mirror1_Voltage = 1000; // Mirror 1 Potential [Volts]
	double Mirror2_Voltage = 2000; // Mirror 2 Potencial [Volts]
	double Field_free1_d = 585e-3; // Grounded tube 1 distance [m]
	double Retard_region1_d = 50e-3; //Mirror 1 distance
	double Retard_region2_d = 97e-3; //Mirror 2 distance
	double Field_free2_d = 311e-3; // Grounded tube 1 distance [m]
	double d0 = 0; // Sample Thickness
	int EXTRACTION_SIZE = 30000;
	double t_Cicle = 1e5;

//------------------------------------------------------------------//
//----------------------------VARIABLES-----------------------------//
//------------------------------------------------------------------//

	int i, j, k, Particle_Index, t_Discrete;
	double t_Dependent_Voltage[EXTRACTION_SIZE], a, b, Eletric_Field, Field_Acceleration, Init_Speed;
	double d, t, v, Eletric_Field_Mirror, Field_Acceleration_Mirror1, v_test, t_Mirror;
	double t_Step_Auxiliar, Field_Acceleration_Mirror2, Eletric_Field_Mirror1, Eletric_Field_Mirror2;
	double v_test2, v_test1;

	//Get 'extraction-pulse2.txt':
	for (i=0;i<EXTRACTION_SIZE;i++)	scanf("%lf %lf %lf", &a, &t_Dependent_Voltage[i], &b);

	Eletric_Field = Extraction_Voltage/Colimator_d;
	Field_Acceleration = (e_Charge*Eletric_Field)/Secondary_I_Mass;
	Eletric_Field_Mirror1 = Mirror1_Voltage/Retard_region1_d;
	Eletric_Field_Mirror2 = Mirror2_Voltage/Retard_region2_d;
	Field_Acceleration_Mirror1 = e_Charge*Eletric_Field_Mirror1/Secondary_I_Mass;
	Field_Acceleration_Mirror2 = e_Charge*Eletric_Field_Mirror2/Secondary_I_Mass;

	for (Particle_Index=0;Particle_Index<1e6;Particle_Index+=10) {
		Init_Speed = gaussian(Central_Velocity);

		if (Particle_Index < t_Init_Pulse) { // Without Extraction Voltage
			d = Init_Speed*(t_Init_Pulse - Particle_Index)*t_Step + d0;
			t = (t_Init_Pulse-Particle_Index)*t_Step;
			if(d>=Colimator_d) {
				t = (Colimator_d - d0)/Init_Speed;
				d = Colimator_d;
			}
		}
		if(d<Colimator_d) {
			v = Init_Speed;
			//t_Step_Auxiliar = t_Step/10;
			while(d<Colimator_d) {
			t_Discrete = Particle_Index + t/t_Step;		
				if(t_Discrete<1e6) {
					t_Discrete = Particle_Index + t/t_Step;
					//printf("%d %lf %d\n", Particle_Index, t/t_Step, t_Discrete);
					t_Discrete = (t_Discrete - t_Init_Pulse);
					
					if(t_Discrete>=30000) t_Discrete = 29999;
					//printf("%d\n", t_Discrete);
					v += t_Dependent_Voltage[t_Discrete]*Field_Acceleration*t_Step;
					d += v*t_Step;
					t += t_Step;
					//printf("%lf \n",t_Step*v,t_Dependent_Voltage[t_Discrete]*Field_Acceleration*t_Step_Auxiliar);
				}	
				if(t_Discrete>=1e6) {
					d += v*t_Step;
					t += t_Step;
				}
			}

			t+= Field_free1_d/v; // Region 2

			v_test1 = sqrt((Mirror1_Voltage)*e_Charge*2/Secondary_I_Mass);
			v_test2 = sqrt((Mirror2_Voltage)*e_Charge*2/Secondary_I_Mass);

			if (v_test1>=v) { // Region 3/4 CASE I
				t_Mirror = 2*v_test1/Field_Acceleration_Mirror1;
				t += t_Mirror;
				t += Field_free2_d/v; // Region 5
				printf("%d %le %lf %lf\n", Particle_Index, t, d, v);
			}

			if (v_test1<=v) { // Region 3/4 CASE II
				t_Mirror = 2*v_test1/Field_Acceleration_Mirror1;
				t += t_Mirror;
				v = v - v_test1;
				if(v_test2>=v) {
					t_Mirror = 2*v_test2/Field_Acceleration_Mirror2;
					t += t_Mirror;
					v = v + v_test1;
					t += Field_free2_d/v; // Region 5
					printf("%d %le %lf %lf\n", Particle_Index, t, d, v);
				}
			}
		}

		if (Particle_Index >= t_Init_Pulse) { // With Extraction Voltage			
			t = 0;
			d = d0; 

			v = Init_Speed;
			//t_Step_Auxiliar = t_Step/10;
			while(d<Colimator_d) {
			t_Discrete = Particle_Index + t/t_Step;		
				if(t_Discrete<1e6) {
					t_Discrete = Particle_Index + t/t_Step;
					//printf("%d %lf %d\n", Particle_Index, t/t_Step, t_Discrete);
					t_Discrete = (t_Discrete - t_Init_Pulse);
					
					if(t_Discrete>=30000) t_Discrete = 29999;
					//printf("%d\n", t_Discrete);
					v += t_Dependent_Voltage[t_Discrete]*Field_Acceleration*t_Step;
					d += v*t_Step;
					t += t_Step;
					//printf("%lf \n",t_Step*v,t_Dependent_Voltage[t_Discrete]*Field_Acceleration*t_Step_Auxiliar);
				}	
				if(t_Discrete>=1e6) {
					d += v*t_Step;
					t += t_Step;
				}
			}

			t+= Field_free1_d/v; // Region 2

			v_test1 = sqrt((Mirror1_Voltage)*e_Charge*2/Secondary_I_Mass);
			v_test2 = sqrt((Mirror2_Voltage)*e_Charge*2/Secondary_I_Mass);

			if (v_test1>=v) { // Region 3/4 CASE I
				t_Mirror = 2*v_test1/Field_Acceleration_Mirror1;
				t += t_Mirror;
				t += Field_free2_d/v; // Region 5
				printf("%d %le %lf %lf\n", Particle_Index, t, d, v);
			}

			if (v_test1<=v) { // Region 3/4 CASE II
				t_Mirror = 2*v_test1/Field_Acceleration_Mirror1;
				t += t_Mirror;
				v = v - v_test1;
				if(v_test2>=v) {
					t_Mirror = 2*v_test2/Field_Acceleration_Mirror2;
					t += t_Mirror;
					v = v + v_test1;
					t += Field_free2_d/v; // Region 5
					printf("%d %le %lf %lf\n", Particle_Index, t, d, v);
				}
			}
		}	

	}

}

double gaussian(double x)
{
	double l=1, PN;
	PN = 1 - 2*(rand()%2);
	return x+0.5*PN*l*log(1.0*rand()/RAND_MAX);
}
