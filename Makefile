OR_TOOLS_TOP=../or-tools
OR_TOOLS_SOURCES=$(OR_TOOLS_TOP)/ortools

include $(OR_TOOLS_TOP)/Makefile

# For debugging uncomment the next line. -isystem prevents warnings rooted in or-tools library appearing in our compilation
CFLAGS := $(CFLAGS) -ggdb -Og -DDEBUG -fsanitize=address -Wall -Wextra -Wshadow -Wunreachable-code -Winit-self -Wmissing-include-dirs -Wswitch-enum -Wfloat-equal -Wundef -isystem$(OR_TOOLS_TOP)/. -isystem$(OR_TOOLS_SOURCES)/gen -isystem$(OR_TOOLS_TOP)/dependencies/install/include -isystem$(OR_TOOLS_TOP)/dependencies/install/include/coin -DUSE_CBC -DUSE_CLP -DUSE_GLOP -DUSE_BOP

.PHONY: all local_clean

all: $(EXE)

%.pb.cc: %.proto
	$(OR_TOOLS_TOP)/dependencies/install/bin/protoc --cpp_out . $<

problem.pb.h: problem.pb.cc

%.o: %.cc %.h
	$(CCC) $(CFLAGS) -c $< -o $@

clustering.o: ./clustering.cc ./problem.pb.h $(OR_TOOLS_SOURCES)/linear_solver/linear_solver.h
	$(CCC) $(CFLAGS) -c ./clustering.cc -o clustering.o

clustering: $(ROUTING_DEPS) clustering.o problem.pb.o $(OR_TOOLS_TOP)/lib/libortools.so
	$(CCC) $(CFLAGS) -g clustering.o problem.pb.o $(OR_TOOLS_LD_FLAGS) \
	-L $(OR_TOOLS_TOP)/lib -lortools -L $(OR_TOOLS_TOP)/dependencies/install/lib -lprotobuf -lglog -lgflags \
	-Wl,-rpath,$(OR_TOOLS_TOP)/dependencies/install/lib -Wl,-rpath,$(OR_TOOLS_TOP)/lib \
	-o clustering

local_clean:
	rm -f *.pb.cc *.pb.h *.o

mrproper: local_clean
	rm -f clustering
