
all: Tree.o Combinatorics.o TreeFit.o
	g++ -g -std=c++11 -o main Tree.o Combinatorics.o TreeFit.o
Tree.o: Node.h Particle.h Tree.h Tree.cpp
	g++ -g -std=c++11 -c Tree.cpp
Combinatorics.o: Combinatorics.h Combinatorics.cpp
	g++ -g -std=c++11 -c Combinatorics.cpp
TreeFit.o: TreeFit.cpp Tree.h Combinatorics.h
	g++ -g -c -std=c++11 TreeFit.cpp
clean:
	rm *.o
	rm main
