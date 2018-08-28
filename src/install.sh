
mkdir -p ./include/R
mkdir -p ./lib
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
echo "loc = R.home('include');write.table(loc,'loc_include.txt',row.names=F,col.names=F,sep=',',quote=F)"  | R --vanilla
locinclude=$(cat loc_include.txt)
cp -r $locinclude/* ./include/R/

mkdir -p ./lib
echo "loc = R.home('lib');write.table(loc,'loc_lib.txt',row.names=F,col.names=F,sep=',',quote=F)"  | R --vanilla
loclib=$(cat loc_lib.txt)
cp -r $loclib/* ./lib/


#Rhome=$(R RHOME)
#cp -r ${Rhome}/include/* ./include/R/
#cp -r ${Rhome}/lib/* ./lib/

ls

#if [[ $OSTYPE == "darwin"* ]]; then

#    install_name_tool -change ~/SFAmix/R/gsl-2.2/lib/libgsl.19.dylib ${PWD}/SFAmix/src/lib/libgsl.19.dylib SFAmix.so
#    install_name_tool -change ~/SFAmix/R/gsl-2.2/lib/libgslcblas.0.dylib ${PWD}/SFAmix/src/lib/libgslcblas.0.dylib SFAmix.so

#fi

#ls


