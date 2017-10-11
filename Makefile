# Makefile

PROGRAM = mTR
OBJS	= main.o handle_one_read.o consensus.o wrap_around_DP.o k_means_clustering.o print_and_feed_one_repeat.o TimeMeasure.o 
CC	= gcc
CFLAGS	= 

.c.o:
	$(CC) $(CFLAGS) -c $<

$(PROGRAM): $(OBJS)
	$(CC) $(OBJS) -lm -o $(PROGRAM)

clean:
	rm $(PROGRAM) $(OBJS)