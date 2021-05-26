# Declaration of variables
CC = mpicxx
EZMATH_DIR = ../ezmath/
FLAGS = -I$(EZMATH_DIR) -DPARALLEL
OPTIMIZATION = -O3
EXEC = conv

EZMATH_SOURCES = Tools.cpp				\
		Point.cpp

SOURCES = main.cpp

ABSOLUTE_EZMATH_SOURCES = $(foreach source,$(EZMATH_SOURCES),$(EZMATH_DIR)$(source))
EZMATH_OBJECTS = $(ABSOLUTE_EZMATH_SOURCES:.cpp=.o)
OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(EZMATH_OBJECTS) $(OBJECTS)
	$(CC) $(EZMATH_OBJECTS) $(OBJECTS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(FLAGS) $(OPTIMIZATION) $< -o $@


# To remove generated files
clean:
	rm -f $(EXEC) ./*.o



