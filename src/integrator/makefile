#makefile for the integrator class and derived variants
include ../make.inc
LIB_PATH=$(LIB_MAIN)

INT_OBJ =  $(LIB_PATH)/rk45.o $(LIB_PATH)/integrator.o $(LIB_PATH)/rk4.o $(LIB_PATH)/euler_sde.o $(LIB_PATH)/example_integrator.o

all: $(INT_OBJ)



$(LIB_PATH)/rk4.o: rk4.cpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/rk4.o
	$(CPP) -c rk4.cpp $(INC_PATH) $(FASTMATH) $(CFLAGS) -o $(LIB_PATH)/rk4.o

$(LIB_PATH)/integrator.o: integrator.cpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/integrator.o
	$(CPP) -c integrator.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/integrator.o

$(LIB_PATH)/rk45.o: rk45.cpp rk45_tmpl.hpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/rk45.o
	$(CPP) -c rk45.cpp $(INC_PATH) $(FASTMATH) $(CFLAGS) -o $(LIB_PATH)/rk45.o

$(LIB_PATH)/lsoda.o: lsoda.cpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/lsoda.o
	$(CPP) -c lsoda.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/lsoda.o 

$(LIB_PATH)/euler_sde.o: euler_sde.cpp euler_sde_tmpl.hpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/euler_sde.o
	$(CPP) -c euler_sde.cpp $(INC_PATH) $(FASTMATH) $(CFLAGS) -o $(LIB_PATH)/euler_sde.o

$(LIB_PATH)/lsoda_tmpl.o: lsoda_tmpl.cpp  $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/lsoda_tmpl.o
	$(CPP) -c lsoda_tmpl.cpp $(INC_PATH) $(FASTMATH) $(CFLAGS) -o $(LIB_PATH)/lsoda_tmpl.o

$(LIB_PATH)/example_integrator.o: example_integrator.cpp example_integrator_tmpl.hpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/example_integrator.o
	$(CPP) -c example_integrator.cpp $(INC_PATH) $(FASTMATH) $(CFLAGS) -o $(LIB_PATH)/example_integrator.o

.PHONY: clean

clean:
	rm -f $(INT_OBJ)

