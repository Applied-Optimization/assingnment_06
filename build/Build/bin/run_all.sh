#!/bin/bash

# Crée le dossier si nécessaire
mkdir -p screen_und_test

echo "Running 10x10 no length..."
./NewtonMethods 0 0 10 10 10000 screen_und_test/results_10x10_no_length.txt

echo "Running 10x10 with length..."
./NewtonMethods 0 1 10 10 10000 screen_und_test/results_10x10_with_length.txt

echo "Running 20x20 no length..."
./NewtonMethods 0 0 20 20 10000 screen_und_test/results_20x20_no_length.txt

echo "Running 20x20 with length..."
./NewtonMethods 0 1 20 20 10000 screen_und_test/results_20x20_with_length.txt

echo "Running 5x5 projected no length..."
./NewtonMethods 1 0 5 5 10000 screen_und_test/results_proj_5x5_no_length.txt

echo "Running 5x5 projected with length..."
./NewtonMethods 1 1 5 5 10000 screen_und_test/results_proj_5x5_with_length.txt

echo "All runs completed. Results stored in 'screen_und_test' folder."

