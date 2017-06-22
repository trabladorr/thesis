all: mult_cycle_test mult_io

mult_cycle_test:
	gcc -o mult_cycle_test mult_cycle_test.c mult_utils.c -O3 -march=native

mult_io:
	gcc -o mult_io mult_io.c mult_utils.c -O3 -march=native

assembly:
	gcc -S mult_cycle_test mult_cycle_test.c mult_utils.c -O3 -march=native
	gcc -S mult_io mult_io.c mult_utils.c -O3 -march=native

clean:
	rm mult_cycle_test &> /dev/null
	rm mult_io &> /dev/null