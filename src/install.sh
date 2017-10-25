

#mkdir -p lib
#mkdir -p include 

#URL='http://www.localmsp.org/gnu/gsl/gsl-2.2.tar.gz'

#wget ${URL}
tar -xzvf gsl-2.2.tar.gz
cd gsl-2.2

./configure --prefix=$PWD
make
make install
cp -r include/gsl ../include/
cp -r lib/* ../lib/

cd ../

#URL='http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz'
#wget ${URL}
#mkdir -p eigen
#tar -xzvf eigen-eigen-67e894c6cd8f.tar.gz
#cp -r eigen-eigen-67e894c6cd8f/Eigen $PWD/include/ 

mkdir -p ./include/R
Rhome=$(R RHOME)
cp -r ${Rhome}/include/* ./include/R/
cp -r ${Rhome}/lib/* ./lib/

ls

#if [[ $OSTYPE == "darwin"* ]]; then

#    install_name_tool -change ~/SFAmix/R/gsl-2.2/lib/libgsl.19.dylib ${PWD}/SFAmix/src/lib/libgsl.19.dylib SFAmix.so
#    install_name_tool -change ~/SFAmix/R/gsl-2.2/lib/libgslcblas.0.dylib ${PWD}/SFAmix/src/lib/libgslcblas.0.dylib SFAmix.so

#fi

#ls


