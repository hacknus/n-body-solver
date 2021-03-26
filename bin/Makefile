############################ -*- Mode: Makefile -*- ###########################
## Makefile ---
###############################################################################


#------------------
CC = CC

# CC compilation options
#-----------------------------
CFLAGS = -I. -I../inc -std=c++17 -O3 -ffast-math

# Linker options
#---------------
LDFLAGS = -O3


# Sources directories
#--------------------
SRCDIR = ../src

# Objects
#------------
OBJ = body.o io.o math_utils.o main.o

main:  $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o main

clean:
	rm -f *.o *.mod

cleanall:
	rm -f *.o *.mod main main_out.*

%.o:    $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $<



