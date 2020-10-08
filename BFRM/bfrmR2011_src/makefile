CXX = g++
#CXX = /usr/local/Cellar/gcc/9.1.0/bin/g++-9
#CXXFLAGS = -O3 -Wno-deprecated -DBFRM_LINUX
CXXFLAGS = -O3 -Wno-deprecated -DBFRM_LINUX


%.o:           	%.cpp
		$(CXX) $(CXXFLAGS) -c $*.cpp 

bfrm_obj = mdp.o bfrmevo.o newran.o extreal.o  Model.o util.o stat.o bfrmresult.o bfrmdata.o bfrmrawdata.o bfrmprior.o newmat1.o newmat2.o newmat3.o newmat4.o newmat5.o newmat6.o newmat7.o newmat8.o newmat9.o newmatex.o bandmat.o submat.o myexcept.o cholesky.o evalue.o hholder.o jacobi.o sort.o svd.o newmatrm.o bfrm.o

bfrm:       	$(bfrm_obj)
		$(CXX) $(CXXFLAGS) -o $@ $(bfrm_obj) -L. -lm
		#$(CXX) $(CXXFLAGS) -static  -o $@ $(bfrm_obj) -L. -lm

newran.o:      	newran.cpp include.h newran.h boolean.h myexcept.h extreal.h

myexcept.o:    	myexcept.cpp include.h boolean.h myexcept.h

extreal.o:     	extreal.cpp include.h boolean.h extreal.h

util.o:		util.cpp stat.cpp util.h stdafx.h newran.h

Model.o:	Model.cpp Model.h

bfrm.o:		bfrm.cpp bfrm.h newran.h newran.cpp tnt.h  stdafx.h Model.h

bfrmevo.o:	bfrmevo.cpp bfrmevo.h bfrm.cpp bfrm.h newran.h newran.cpp

util.o:		util.h util.cpp newran.cpp newran.h

stat.o:		util.h stat.cpp

bfrmresult.o:	bfrmresult.h bfrmresult.cpp

bfrmrawdata.o:	bfrmdata.h bfrmdata.cpp bfrmrawdata.h bfrmrawdata.cpp

bfrmdata.o:	bfrmdata.h bfrmdata.cpp

bfrmprior.o:	bfrmprior.h bfrmprior.cpp newran.cpp newran.h Model.h Model.cpp stdafx.h util.h util.cpp

mdp.o:		mdp.cpp mdp.h newran.h stdafx.h

newmat1.o:     	newmat1.cpp newmat.h nminclude.h boolean.h myexcept.h

newmat2.o:     	newmat2.cpp nminclude.h newmat.h newmatrc.h boolean.h myexcept.h controlw.h

newmat3.o:     	newmat3.cpp nminclude.h newmat.h newmatrc.h boolean.h myexcept.h controlw.h

newmat4.o:     	newmat4.cpp nminclude.h newmat.h newmatrc.h boolean.h myexcept.h controlw.h

newmat5.o:     	newmat5.cpp nminclude.h newmat.h newmatrc.h boolean.h myexcept.h controlw.h

newmat6.o:     	newmat6.cpp nminclude.h newmat.h newmatrc.h boolean.h myexcept.h controlw.h

newmat7.o:     	newmat7.cpp nminclude.h newmat.h newmatrc.h boolean.h myexcept.h controlw.h

newmat8.o:     	newmat8.cpp nminclude.h newmat.h newmatrc.h precisio.h boolean.h myexcept.h controlw.h

newmatex.o:    	newmatex.cpp nminclude.h newmat.h boolean.h myexcept.h

bandmat.o:     	bandmat.cpp nminclude.h newmat.h newmatrc.h boolean.h myexcept.h controlw.h

submat.o:      	submat.cpp nminclude.h newmat.h newmatrc.h boolean.h myexcept.h controlw.h

cholesky.o:    	cholesky.cpp nminclude.h newmat.h boolean.h myexcept.h

evalue.o:      	evalue.cpp nminclude.h newmatap.h newmatrm.h precisio.h newmat.h boolean.h myexcept.h

fft.o:         	fft.cpp nminclude.h newmatap.h newmat.h boolean.h myexcept.h

hholder.o:     	hholder.cpp nminclude.h newmatap.h newmat.h boolean.h myexcept.h

jacobi.o:      	jacobi.cpp nminclude.h newmatap.h precisio.h newmatrm.h newmat.h boolean.h myexcept.h

newfft.o:      	newfft.cpp nmnewmatap.h newmat.h include.h boolean.h myexcept.h

sort.o:        	sort.cpp nminclude.h newmatap.h newmat.h boolean.h myexcept.h

svd.o:         	svd.cpp nminclude.h newmatap.h newmatrm.h precisio.h newmat.h boolean.h myexcept.h

newmatrm.o:    	newmatrm.cpp newmat.h newmatrm.h nminclude.h boolean.h myexcept.h

newmat9.o:     	newmat9.cpp nminclude.h newmat.h newmatio.h newmatrc.h boolean.h myexcept.h controlw.h



