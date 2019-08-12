/**
 * High-performance Computing
 *
 * Example code demonstrating the use of the JPEG library to read and write
 * image files to/from arrays.
 */
#include <iostream>
#include <iostream>
using namespace std;

#include <boost/gil/extension/io/jpeg_io.hpp>
using namespace boost::gil;

int main() {
    // Read the image
    rgb8_image_t img;
    jpeg_read_image("plane.jpg", img);

    const int n = img.width();
    const int m = img.height();

    cout << "Image dimensions: " << n << " x " << m << endl;

    // Extract the pixel data to three matrices
    double* r = new double[n*m];
    double* g = new double[n*m];
    double* b = new double[n*m];
    auto imgv = const_view(img);
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            r[x*m+y] = (int)at_c<0>(imgv(x,y));
            g[x*m+y] = (int)at_c<1>(imgv(x,y));
            b[x*m+y] = (int)at_c<2>(imgv(x,y));
        }
    }

    // ... do stuff here ...

    // Generate a new image to write out
    rgb8_image_t out(n, m);
    auto outv = view(out);
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            outv(x,y) = rgb8_pixel_t(r[x*m+y], g[x*m+y], b[x*m+y]);
        }
    }
    jpeg_write_view("output.jpg", outv);

    delete[] r;
    delete[] g;
    delete[] b;
}
