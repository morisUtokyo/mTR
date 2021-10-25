# Makefile

PROGRAM = mTR
OBJS	= main.o handle_one_read.o consensus.o wrap_around_DP.o handle_one_file.o chaining.o fill_directional_index.o  
CC	= mpicc
CPP	= mpic++
CFLAGS	= -std=c99 -fPIC -fcommon

.cpp.o:
	$(CPP) -c $< 
.c.o:
	$(CC) $(CFLAGS) -c $<

# g++ must be used to link libraries required
$(PROGRAM): $(OBJS)
	$(CPP) $(OBJS) -o $(PROGRAM)
clean:
	rm $(PROGRAM) $(OBJS)
