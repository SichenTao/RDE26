CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wno-unused-parameter
LDLIBS ?=

TARGET := rde2026_mop_r05_t0070
SRC := rde2026_mop_r05_t0070.cpp

.PHONY: all clean smoke

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

smoke: $(TARGET)
	./$(TARGET) --problem 1 --run 1 --N 30 --M 3 --D 7 --maxFE 2000 --save 21 --seed 2025171001 --pf-root ../../../benchmark/SEC2018_MaOP_M3_D10/POF --out-dir ../results/smoke

clean:
	rm -f $(TARGET)
