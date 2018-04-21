
all: Tree.o TreeFit.o
	g++ -g -std=c++11 -o main Tree.o TreeFit.o
Tree.o: Node.h Particle.h Tree.h Tree.cpp
	g++ -g -std=c++11 -c Tree.cpp
TreeFit.o: TreeFit.cpp
	g++ -g -c -std=c++11 TreeFit.cpp
clean:
	rm *.o
	rm main
