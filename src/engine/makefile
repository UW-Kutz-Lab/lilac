include ../make.inc
LIB_PATH=$(LIB_MAIN)
ENG_OBJ= $(LIB_PATH)/engine_main.o $(LIB_PATH)/item.o $(LIB_PATH)/engine_run.o $(LIB_PATH)/engine_write.o $(LIB_PATH)/variable.o  $(LIB_PATH)/item_dim.o $(LIB_PATH)/input.o $(LIB_PATH)/item_wrapper.o $(LIB_PATH)/async_write.o
all: $(ENG_OBJ)

$(LIB_PATH)/engine_main.o:engine.cpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/engine_main.o
	$(CPP) -c engine.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/engine_main.o

$(LIB_PATH)/item.o:item.cpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/item.o
	$(CPP) -c item.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/item.o

$(LIB_PATH)/engine_run.o: engine_run.cpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/engine_run.o
	$(CPP) -c engine_run.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/engine_run.o

$(LIB_PATH)/engine_write.o:engine_write.cpp ../writer/data.hpp $(INLINE_HEADERS) 
	@rm -f $(LIB_PATH)/engine_write.o
	$(CPP) -c engine_write.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/engine_write.o

$(LIB_PATH)/variable.o:variable.cpp $(INLINE_HEADERS) 
	@rm -f $(LIB_PATH)/variable.o
	$(CPP) -c variable.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/variable.o

$(LIB_PATH)/item_dim.o:item_dim.cpp $(INLINE_HEADERS)
	@rm -f $(LIB_PATH)/item_dim.o
	$(CPP) -c item_dim.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/item_dim.o

$(LIB_PATH)/input.o:input.cpp $(INLINE_HEADERS) 
	@rm -f $(LIB_PATH)/input.o
	$(CPP) -c input.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/input.o

$(LIB_PATH)/item_wrapper.o:item_wrapper.cpp $(INLINE_HEADERS) 
	@rm -f $(LIB_PATH)/item_wrapper.o
	$(CPP) -c item_wrapper.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/item_wrapper.o


$(LIB_PATH)/list.o:list.cpp $(INLINE_HEADERS) 
	@rm -f $(LIB_PATH)/list.o
	$(CPP) -c list.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/list.o

$(LIB_PATH)/async_write.o:async_write.cpp $(INLINE_HEADERS) 
	@rm -f $(LIB_PATH)/async_write.o
	$(CPP) -c async_write.cpp $(INC_PATH) $(CFLAGS) -o $(LIB_PATH)/async_write.o


.PHONY: clean

clean:
	rm -f $(ENG_OBJ)
