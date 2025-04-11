CXX = g++
CXXFLAGS = -O2 -std=c++11
TARGET = test_sequential
SRCS = sequential.cpp test_sequential.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

clean:
	rm -f $(TARGET)