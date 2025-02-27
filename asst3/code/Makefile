APP_NAME=wireroute

OBJS=wireroute.o validate.o

CXX = g++
CXXFLAGS = -Wall -O3 -std=c++17 -m64 -I. -fopenmp -Wno-unknown-pragmas -pg

all: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class
	/bin/rm inputs/testinput/*_occupancy_*.txt
	/bin/rm inputs/testinput/*_wires_*.txt
	/bin/rm inputs/timeinput/*_occupancy_*.txt
	/bin/rm inputs/timeinput/*_wires_*.txt
