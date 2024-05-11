all: main.exe

main.exe: main.o kDTree.o Dataset.o
	gcc -o main main.cpp kDTree.cpp Dataset.o -I . -std=c++11

clean:
	del  main