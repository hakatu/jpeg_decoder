rm -r build
mkdir build
cd build
cmake ..
make
./jpeg_decode ../../test/space.jpg bitmap.bmp
