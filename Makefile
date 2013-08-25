CC = g++
CFLAGS  = -Wall
LDFLAGS = -lboost_program_options

all: drstats

drstats: dvect.o stats.o main.o
	$(CC) $(CFLAGS) $(LDFLAGS) dvect.o stats.o main.o -o drstats

main.o: stats.h main.cpp
	$(CC) $(CFLAGS) -c main.cpp

stats.o: dvect.h stats.h stats.cpp
	$(CC) $(CFLAFS) -c stats.cpp

dvect.o: dvect.h dvect.cpp
	$(CC) $(CFLAGS) -c dvect.cpp

clean:
	rm drstats *.o *.gch
