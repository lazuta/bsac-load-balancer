all: run
run: simulator.o graph.o logger.o
	g++ -orun.exe simulator.o graph.o logger.o run.cpp
simulator.o: graph.o logger.o simulator.cpp
	g++ -c -osimulator.o simulator.cpp
graph.o: preflow_graph.cpp
	g++ -c -ograph.o preflow_graph.cpp
logger.o: logger.cpp
	g++ -c -ologger.o logger.cpp

