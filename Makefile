
CXXFLAGS += -std=c++20 -pedantic
CXXFLAGS += -Wall -Wextra -Werror

main: CXXFLAGS += -O3 -march=native -DNDEBUG
debug: CXXFLAGS += -O2 -g

main: *.cpp *.hpp params.hpp
	$(CXX) $(CXXFLAGS) *.cpp -o $@

debug: *.cpp *.hpp params.hpp
	$(CXX) $(CXXFLAGS) *.cpp -o $@

SAGE ?= sage
params.hpp: params.sage
	$(SAGE) $< > $@

.PHONY: clean
clean:
	rm -f main debug params.hpp
	rm -f *.sage.py

