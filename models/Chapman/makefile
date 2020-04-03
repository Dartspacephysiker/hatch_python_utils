# Makefile for building Tsyganeko and compliling it with f2py

# Declare default compilers
LOCAL_F2PY = f2py
LOCAL_F77  = gfortran
# Use default compilers unless specified by command line input
#echo 'Use F2PY and F77 to specify the fortran-to-python and fortran complilers'
F2PY = $(LOCAL_F2PY)
F77  = $(LOCAL_F77)

# Declare optional and required compilation flags
# OPT_FLAGS = -O2 -fbacktrace -fno-automatic -fPIC
REQ_FLAGS = -w

# Declare source code
SOURCES = chapman.for
# SRC = geopack.for
OBJ = $(SOURCES:.for=.o) 

PYF = $(SRC:.for=.pyf)
PYO = $(PYF:.pyf=.so)

# Define operations
all:
        # $(F77) $(REQ_FLAGS) $(OPT_FLAGS) -c $(SRC) $(SOURCES)
	$(F2PY) --f77flags="$(REQ_FLAGS)" -c -m chapman $(PYF) $(SOURCES)

clean:
	rm -f $(OBJ) $(PYO) chapman.cpython*.so

# Declare test operations
# test:
# 	$(F2PY) $(SRC) -m tsygFort -h geopack_test.pyf $(SOURCES)

# test_clean:
# 	rm -f geopack_test.pyf
