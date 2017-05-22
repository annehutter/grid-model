SOURCES := 	./src/main.c 
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := create_parameter_file

include ../common.mk

# OPTIMIZE = -O3 -ftree-vectorize
# WARNING = -Wall -Wextra -Wshadow -g


.PHONY: all clean

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
