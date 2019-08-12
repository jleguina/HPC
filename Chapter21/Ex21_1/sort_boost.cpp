/* 
 * High-performance Computing
 *
 * Insertion sort with boost program options.
 */
#include <iostream>
#include <cstdlib>
using namespace std;

#include <boost/program_options.hpp>

// An alias to reduce typing
namespace po = boost::program_options;

// Insert a value into the correct position in array
void insert(float a[], int p, float val, bool desc) {
	if (desc) {
        // Move items up array until the right place is reached
        while (p > 0 && val > a[p-1]) {
            a[p] = a[p-1];
            p--;
        }
	}
    else {
        // Move items up array until the right place is reached
        while (p > 0 && val < a[p-1]) {
            a[p] = a[p-1];
            p--;
        }
    }
    // Put 'val' in the right place
	a[p] = val;
}

// Perform an insertion sort on an array
void insertionSort(float a[], int length, bool desc) {
	for (int i = 0; i < length; i++) {
		insert(a,i,a[i],desc);
	}
}

int main(int argc, char* argv[]) {
    po::options_description opts(
		"Sorts a list of random numbers using the insertion sort algorithm.");
    opts.add_options()
        ("size", po::value<int>()->default_value(10),
				 "Size of vector to sort.")
        ("min",  po::value<int>()->default_value(0),
				 "Minimum value of numbers to generate.")
        ("max",  po::value<int>()->default_value(10),
				 "Maximum value of numbers to generate.")
		("descending", "Indicate the array should be reversed.")
        ("help",       "Print help message.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << "Performs an insertion sort algorithm on an array of "
             << "random numbers." << endl;
        cout << opts << endl;
        return 0;
    }

    const int n = vm["size"].as<int>();
    const int range = vm["max"].as<int>() - vm["min"].as<int>();
	const int offset = vm["min"].as<int>();
	const bool desc = vm.count("descending");

    float *a = new float[n];

    srand(time(0));

    // Generate and print random numbers
    cout << "Generated random numbers: " << endl;
    for (int i = 0; i < n; i++) {
		a[i] = float(rand()) / RAND_MAX * range + offset;
		cout << a[i] << endl;
    }
	cout << endl;

	// Sort array
	insertionSort(a,n,desc);

	// Print sorted array
    cout << "Sorted random numbers: " << endl;
	for (int i = 0; i < n; i++) {
		cout << a[i] << endl;
	}

    delete[] a;

    return 0;
}

