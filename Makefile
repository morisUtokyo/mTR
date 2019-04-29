# Makefile

PROGRAM = mTR
OBJS	= main.o handle_one_read.o consensus.o wrap_around_DP.o print_and_feed_one_repeat.o handle_one_file.o chaining.o  
CC	= gcc
CPP	= g++
CFLAGS	= -std=c99

.cpp.o:
	$(CPP) -c $< 
.c.o:
	$(CC) $(CFLAGS) -c $<

# g++ must be used to link libraries required
$(PROGRAM): $(OBJS)
	$(CPP) $(OBJS) -o $(PROGRAM)
clean:
	rm $(PROGRAM) $(OBJS)
