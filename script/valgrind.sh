#!/bin/sh

# use valgrind memcheck
valgrind --tool=memcheck --leak-check=full --track-origins=yes $*
