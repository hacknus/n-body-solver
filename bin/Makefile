CC = CC
CFLAGS = -I. -I../inc -O3 -ffast-math
SRCDIR = ../src

all: main.o body.o math_utils.o io.o
	$(CC) $(CFLAGS) -o main main.o body.o math_utils.o io.o

main.o: $(SRCDIR)/main.cpp
	$(CC) $(CFLAGS) -c $(SRCDIR)/main.cpp

body.o: $(SRCDIR)/body.cpp
	$(CC) $(CFLAGS) -c $(SRCDIR)/body.cpp

io.o: $(SRCDIR)/io.cpp
	$(CC) $(CFLAGS) -c $(SRCDIR)/io.cpp

math_utils.o: $(SRCDIR)/math_utils.cpp
	$(CC) $(CFLAGS) -c $(SRCDIR)/math_utils.cpp

clean:
	rm main main.o body.o math_utils.o io.o