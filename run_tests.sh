#!/bin/bash

# Loop over all .txt files in the tests/unit_scenes directory
for file in ./tests/unit_scenes/*.txt; do
    echo "Running test: $file"
    ./bin/main.exe "$file"
    echo ""  # Print a newline for readability
done
