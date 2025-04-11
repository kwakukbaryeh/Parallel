CXX = g++
CXXFLAGS = -O2 -std=c++11
TARGET = test_large_sequential
SRCS = sequential.cpp test_large_sequential.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

clean:
	rm -f $(TARGET)
