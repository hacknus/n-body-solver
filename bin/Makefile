############################ -*- Mode: Makefile -*- ###########################
## Makefile ---
###############################################################################


#------------------
# CC=CC for cray system on CSCS
# CC=mpic++ on macOS
CC = mpic++

# CC compilation options
#-----------------------------
CFLAGS = -I../inc -O3 -std=c++17 -ffast-math

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
	$(CC) $(LDFLAGS) $(OBJ) -o nbody

clean:
	rm -f *.o *.mod

cleanall:
	rm -f *.o *.mod nbody nbody.out nbody.err

%.o:    $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $<



