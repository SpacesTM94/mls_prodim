
mls.out: points.cpp
	g++ -std=c++17 points.cpp -o $@ -lstdc++

mls_test.out: points.cpp
	g++ -std=c++17 points.cpp -o $@ -DRUN_TESTS -lstdc++

clear:
	rm mls.out
	rm mls_test.out

all: mls.out mls_test.out
