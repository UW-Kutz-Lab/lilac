#makefile for mpi components

include ../make.inc

LIB_PATH=$(LIB_MAIN)

MPI_OBJ = $(LIB_PATH)/mpi_controller.o

all: $(MPI_OBJ)

$(LIB_PATH)/mpi_controller.o: mpi_controller.cpp mpi_controller.h
	@rm -f $(LIB_PATH)/mpi_controller.o
	$(CPP) -c mpi_controller.cpp $(INC_PATH)  -o $(LIB_PATH)/mpi_controller.o  $(CFLAGS) $(FASTMATH)
 
