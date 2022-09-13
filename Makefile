CC=g++
CFLAGS=-O3 -std=c++17
DEPS = 
OBJ = collatz.o 
LIBS = -lboost_chrono

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

collatz: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
