CC=g++
CFLAGS=-O3
DEPS = 
OBJ = collatz.o 

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

collatz: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
