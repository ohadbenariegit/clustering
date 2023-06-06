#!/bin/bash
# Script to compile and execute a c program
gcc -ansi -Wall -Wextra -Werror -pedantic-errors common/helper_methods.c kmeans/kmeans.c spkmeans.c -lm -o spkmeans