# Makefile

PROGRAM = mTR
OBJS	= main.o handle_one_file.o handle_one_read.o consensus.o wrap_around_DP.o print_and_feed_one_repeat.o 
CC	= gcc
CFLAGS	= -g -std=c99

.c.o:
	$(CC) $(CFLAGS) -c $<

$(PROGRAM): $(OBJS)
	$(CC) $(OBJS) -o $(PROGRAM)

clean:
	rm $(PROGRAM) $(OBJS)
