#!/bin/sh
function timeout() { perl -e 'alarm shift; exec @ARGV' "$@"; }
timeout 2 ./tcp2run -n 2 -s 1 -r L1 ../../examples/resize 1 | LC_ALL='C' sort > test-resize-2-s1-r1.out
cmp test-resize-2-s1-r1.out "$(dirname -- "${0}")/test-resize-2-s1-r1.expected"