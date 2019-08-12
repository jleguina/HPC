#include <iostream>
#define G -9.81
#include <math.h>

 using namespace std;

 int main()
 {

float vWeight;
float vDragCoeff;

 cout << "Enter mass (in kg): ";
 cin >> vWeight ;

 while (vWeight <0) {
 cout << "Mass must be a positive number. Try again: ";
 cin >> vWeight ;
 }

  cout << "Enter drag coefficient (in kg/s): ";
 cin >> vDragCoeff ;

 float vTime1=3;
 float vVelocity1= (G*vWeight/vDragCoeff)*(pow(2.41,-vDragCoeff*vTime1/vWeight)-1);

 float vTime2=5;
 float vVelocity2= (G*vWeight/vDragCoeff)*(pow(2.41,-vDragCoeff*vTime2/vWeight)-1);

 float vTime3=10;
 float vVelocity3= (G*vWeight/vDragCoeff)*(pow(2.41,-vDragCoeff*vTime3/vWeight)-1);


 cout << "The velocity at t=3 is: " << vVelocity1 << "m/s" << endl;
 cout << "The velocity at t=5 is: " << vVelocity2 << "m/s" << endl;
 cout << "The velocity at t=10 is: " << vVelocity3 << "m/s" << endl;



 }
