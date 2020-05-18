# compilor
CC=g++ -std=c++11
# mode
ifndef MODE
	MODE=RELEASE
endif

CFLAGS=-Wall
ifeq ($(MODE),DEBUG)
CFLAGS=-Wall -g
endif

objects=main.o tools.o debug.o

AVE:$(objects)
	$(CC) $(CFLAGS) $(objects) -o ave

main.o:main.cpp
	$(CC) $(CFLAGS) -c main.cpp
tools.o:tools.cpp
	$(CC) $(CFLAGS) -c tools.cpp
#geometry.o:geometry.cpp
#	$(CC) $(CFLAGS) geometry.cpp
#ellipse.o:ellipse.cpp
#	$(CC) $(CFLAGS) ellipse.cpp
debug.o:debug.cpp
	$(CC) $(CFLAGS) -c debug.cpp
clear:
	rm ./*.o
