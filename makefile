CC=g++
CFLAGS=  -Wall -Wextra  -Ofast -std=c++11  -pthread -pipe -fopenmp
LDFLAGS=-pthread -fopenmp


EXEC=btt

all: $(EXEC)

btt:   btt.o
	$(CC) -o $@ $^ $(LDFLAGS)

btt.o:   btt.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
