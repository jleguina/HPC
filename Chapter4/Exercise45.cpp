#include <iostream>
#define G -9.81
 using namespace std;

 int main()
 {

float vHeight;
float vTimeInitial;
float vTimeFinal;

 cout << "Enter initial height: ";
 cin >> vHeight ;

 cout << "Enter initial time: ";
 cin >> vTimeInitial ;

 cout << "Enter final time: ";
 cin >> vTimeFinal ;

 float vVelocity= (G)*(vTimeFinal-vTimeInitial);

 cout << "The velocity is: " << vVelocity << "m/s" << endl;



 }
