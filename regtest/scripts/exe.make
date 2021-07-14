
exe:
	$(CXX) -c $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) *.cpp
	$(LD) *.o -o $@ $(PLUMED_LOAD)

exe-fortran:
	$(FC) -c $(PLUMED_FORTRAN) *.f90
	$(FC) *.o -o exe $(PLUMED_LOAD)

print-fortran:
	@echo FC=$(FC)
