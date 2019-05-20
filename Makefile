GPU_ARCH ?= sm_13

example_dd: example_dd.cu dbldbl.h
	nvcc -I . -o example_dd -arch=$(GPU_ARCH) example_dd.cu

clean:
	rm -f example_dd *~
