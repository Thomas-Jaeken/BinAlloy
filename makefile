CC = gcc
CFLAGS = -O3 -Wall -fsanitize=address -fno-omit-frame-pointer -g
LIBS = -lm -O3 -lgsl -lgslcblas

HEADERS = 
OBJECTS = main.o
PROGRAM = mean_field
OBJECTS2 = MCmain.o
PROGRAM2 = monte_carlo_s


%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	rm -f *.o
	touch *.c

mc: $(PROGRAM2)

$(PROGRAM2): $(OBJECTS2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	rm -f *.o
	touch *.c