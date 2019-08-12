#include <iostream>

 using namespace std;

 int main()
 {

float vRadius;

 cout << "Enter radius of sphere: ";
 cin >> vRadius ;

 float vVolume= 4*3.1415*vRadius*vRadius*vRadius/3;
 float vSurface=4*3.1415*vRadius*vRadius;

 cout << "Volume of sphere is: " << vVolume << endl;
 cout << "Surface of sphere is: " << vSurface << endl;



 }
