#include <iostream>
#include <cmath>

using namespace std;

int main()
{

unsigned int n;
unsigned int k;
int fact=1 ;

cout << "Input n (>0 integer):" << endl;
cin >> n ;


for (k=1; k<=n; k++){
fact *= k ;

}

cout << fact << endl;
}
