CXX=mpic++ 
#CXX = g++
#CXX = tau_cxx.sh

APP = hedonics-v2.10
#APP = hedonics-runs

hedonics: $(APP).o hedonics-utilities-v2.10.o 
	$(CXX) -o hedonics $(APP).o hedonics-utilities-v2.10.o

$(APP).o: $(APP).cpp
	$(CXX) -c $(APP).cpp

hedonics-utilities-v2.10.o: hedonics-utilities-v2.10.cpp
	$(CXX) -c hedonics-utilities-v2.10.cpp

clean:
	rm hedonics $(APP).o hedonics-utilities-v2.10.o
