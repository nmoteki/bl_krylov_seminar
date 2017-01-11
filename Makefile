TARGET = bl_krylov_demo
CXX = g++-6
CXXFLAGS = -O3
CPPFLAGS = -DEIGEN_NO_DEBUG -DNDEBUG
INCLUDES = -I/Users/moteki/eigen_3_2_10
$(TARGET) : bl_krylov_demo.o bl_bicg.o bl_bicg_rq.o qr_reduced.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -o $@ bl_krylov_demo.o bl_bicg.o bl_bicg_rq.o qr_reduced.o

bl_krylov_demo.o : bl_krylov_demo.cpp bl_krylov_solvers.hpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_krylov_demo.cpp

bl_bicg.o : bl_bicg.cpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_bicg.cpp

bl_bicg_rq.o : bl_bicg_rq.cpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_bicg_rq.cpp

qr_reduced.o: qr_reduced.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c qr_reduced.cpp

clean:
	rm -f *.o $(TARGET)
