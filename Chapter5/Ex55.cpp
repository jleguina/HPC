#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>


using namespace std;

int main()
{
	
ofstream vOut("sine.txt", ios::out| ios::trunc);
vOut.precision(5);

int i;
double s [100];

for (i=0; i<100; i++) {
	
s[i]=sin(2.*M_PI*i/100.);

vOut << setw(10) << s[i] << setw(10) << i << endl;

cout << s[i] << endl;

}

vOut.close();



}
