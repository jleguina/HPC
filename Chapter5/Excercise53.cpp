#include <iostream>
#include <cmath>

using namespace std;

int main()
{

int k;
float pi2;
float pi4;

for (k=0; k<100; k++){
pi4 = pow(-1,k)/(2*k+1);

}

for (k=1; k<100; k++){
pi2 = 4*k*k/(4*k*k-1);

}
cout << 4*pi4 << endl;
cout << 2*pi2 << endl;
}
