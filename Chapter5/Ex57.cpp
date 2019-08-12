#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>


using namespace std;

int main()
{
	



float g=9.81;
float c=16;
float m=65;

float i=1;
float dt=0.1;
float v=0;
float t=0;

ofstream vOut("velocity.txt", ios::out| ios::trunc);
vOut.precision(5);

for (i=2;i<100000;i++){
	t=dt*i;
	v=v+(g-(c*v/m))*(dt);
	
vOut << setw(10) << t << setw(10) << v << endl;

cout << v << endl;
	//cout << t << endl;
}
	




vOut.close();
	
	

}
