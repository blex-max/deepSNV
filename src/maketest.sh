#!/usr/bin/env bash

# use bear -- before make (not clean) to get compile_commands.json
make -f Makefile.test clean && make -f Makefile.test && ./build/pileup-test-exec
