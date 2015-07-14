CXX := g++
#CXX := clang++
#CXX := icpc #-gxx-name=/usr/bin/g++-4.7

CXXFLAGS := -std=c++11 -O2 -fpic -I../../src -I. -g -Wall -Wextra -Wconversion
LDFLAGS := -shared -g


SOURCES := linear.cc fixed_power.cc templated_power.cc variable_power.cc sinc.cc
OBJECTS  := $(patsubst %.cc,%.o,$(SOURCES))
SHARED_OBJECTS := $(patsubst %.cc,%.so,$(SOURCES))


.PHONY: all
all: $(SHARED_OBJECTS)


.PHONY: clean
clean:
	$(RM) *.so *.o


.PHONY: distclean
distclean: clean
	$(RM) *~


%.so: %.o
	$(CXX) $(LDFLAGS) -Wl,-soname,$@ -o $@ $^


%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $<


# All-in-one library
userweights.so: $(OBJECTS)


sinc.o: sinc.cc root_finding.hh


# User-weight `tabular' is not part of SOURCES, because it relies on a
# GSL header, which may not available on every system.
tabular.so: tabular.o interpolator.o
tabular.o: tabular.cc interpolator.hh

interpolator.cc: interpolator.hh
