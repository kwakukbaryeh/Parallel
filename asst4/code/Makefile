APP_NAME=wireroute

OBJS=wireroute.o validate.o

CXX = mpic++
CXXFLAGS = -Wall -Wextra -O3 -std=c++20 -I.

all: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class
