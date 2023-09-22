CC=g++
CFLAGS=-O3 -std=c++17 -Wall
DEPS = 
OBJ = collatz_gmp.o 
LIBS = -lboost_chrono -lgmp

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

collatz_gmp: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

