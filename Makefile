CC = g++
CCFLAGS = -O4
all: simulator.o graph.o utils.o
	$(CC) $(CCFLAGS) -orun.exe simulator.o graph.o logger.o rational.o run.cpp
simulator.o: graph.o simulator.cpp
	$(CC) $(CCFLAGS) -c -osimulator.o simulator.cpp
graph.o: preflow_graph.cpp utils.o
	$(CC) $(CCFLAGS) -c -ograph.o preflow_graph.cpp
utils.o: 
	$(CC) $(CCFLAGS) -c -ologger.o utils/logger.cpp
	$(CC) $(CCFLAGS) -c -orational.o utils/rational.cpp
clean:
	$(RM) *.o run.exe
