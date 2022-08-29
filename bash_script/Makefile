# CC = /usr/local/bin/gcc-11
# CFLAGS = -O3 -std=c++11 -march=native -mtune=native -fno-stack-protector -fomit-frame-pointer
# LDFLAGS = -L /usr/local/opt/gcc/lib/gcc/11 -L /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lSystem -lstdc++

CC = clang++
CFLAGS = -std=c++11 -Wall -Wextra -pedantic -O3 -march=native -mtune=native -flto -fno-stack-protector -fomit-frame-pointer
# CFLAGS = -std=c++11 -Wall -Wextra -pedantic -g -fno-stack-protector -fomit-frame-pointer
LDFLAGS =
# OMPFLAGS =
# OMPFLAGS = -Xclang -fopenmp -lomp

all: CatForce active_region_occurs compute_flipper_filters
CatForce: CatForce.cpp LifeAPI.h
	$(CC) $(CFLAGS) -o CatForce CatForce.cpp $(LDFLAGS)
active_region_occurs: active_region_occurs.cpp LifeAPI.h
	$(CC) $(CFLAGS) -o active_region_occurs active_region_occurs.cpp $(LDFLAGS)
compute_flipper_filters: compute_flipper_filters.cpp LifeAPI.h
	$(CC) $(CFLAGS) -o compute_flipper_filters compute_flipper_filters.cpp $(LDFLAGS)

