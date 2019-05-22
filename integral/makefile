TARGET = int

$TARGET:
	g++ -std=c++11 -I/home/ydzhao/cuth/src/Cuba-4.2 -I/home/ydzhao/cuth/install/boost/include -O3 -c $(TARGET).cc -o $(TARGET).o
	g++ -L/home/ydzhao/cuth/src/Cuba-4.2 -L/home/ydzhao/cuth/install/boost/lib $(TARGET).o -lcuba -o $(TARGET) -lboost_program_options -lgsl

