#!/bin/bash

# Copyright (C) 2018 by Keith Pedersen (Keith.David.Pedersen@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# This script determines GCC vector extension flags on a Linux based OS

flags=""

if grep -q sse4_2 /proc/cpuinfo; then
	flags+=' -msse4.2'
elif grep -q sse4_1 /proc/cpuinfo; then
	flags+=' -msse4.1'
elif grep -q sse4 /proc/cpuinfo; then
	flags+=' -msse4'
elif grep -q sse3 /proc/cpuinfo; then
	flags+=' -msse3'
elif grep -q sse2 /proc/cpuinfo; then
	flags+=' -msse2'
elif grep -q sse /proc/cpuinfo; then
	flags+=' -msse'
fi

if grep -q avx2 /proc/cpuinfo; then
	flags+=' -mavx2'
elif grep -q avx /proc/cpuinfo; then
	flags+=' -mavx'
fi

echo $flags
