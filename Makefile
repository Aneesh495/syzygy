GHC ?= ghc
GHCFLAGS ?= -O2 -Wall -Wcompat -i src -i app -i test -i bench -outputdir dist/build
CC ?= clang
CFLAGS ?= -O3 -fPIC -Wall -Wextra -std=c11

.PHONY: all build c-lib test bench check-dashes clean repl

all: check-dashes c-lib build test

build:
	@mkdir -p bin dist/build
	$(GHC) $(GHCFLAGS) app/Main.hs -o bin/syzygy

c-lib:
	@mkdir -p dist/build
	$(CC) $(CFLAGS) -shared c/syzygy_kernel.c -o dist/build/libsyzygy_c.dylib 2>/dev/null || \
	$(CC) $(CFLAGS) -shared c/syzygy_kernel.c -o dist/build/libsyzygy_c.so 2>/dev/null || true

test:
	@mkdir -p bin dist/build
	$(GHC) $(GHCFLAGS) test/Spec.hs -o bin/test_runner
	./bin/test_runner

bench:
	@mkdir -p bin dist/build
	$(GHC) $(GHCFLAGS) bench/Bench.hs -o bin/bench_runner
	./bin/bench_runner

bench-c:
	@mkdir -p bin
	$(CC) $(CFLAGS) c/bench_kernel.c c/syzygy_kernel.c -o bin/bench_kernel
	./bin/bench_kernel

check-dashes:
	python3 tools/check_no_em_dashes.py

clean:
	rm -rf dist/ bin/
	find . -type f \( -name "*.o" -o -name "*.hi" -o -name "*.dyn_o" -o -name "*.dyn_hi" \) -delete

repl:
	ghci -i src
