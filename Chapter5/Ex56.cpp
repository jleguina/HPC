#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>


using namespace std;

int main()
{
	
ifstream vMyFile("sine.txt");

string vTemp;
float var1, var2;

if (vMyFile.good()) { 			// Check file opened ok
while (true) { 					// Keep trying...
vMyFile >> var1 >> var2; 				// ... to read file
if (vMyFile.eof()) { 			// Test if end of file reached
break; 							// ... and stop if it is
}
cout << var1 << endl; 			// Print what we read
cout << var2 << endl; 
}
vMyFile.close();
}
else {
cout << "Failed to open file" << endl;
}


}
