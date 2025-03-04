#!/bin/bash

# Directory containing the markdown files
DIR="./"

# Replace \textbf{...} with **...**
find "$DIR" -name "*.md" -exec sed -i 's/\\textbf{\([^}]*\)}/**\1**/g' {} +

# Replace \_ with _
find "$DIR" -name "*.md" -exec sed -i 's/\\_/_/g' {} +

# Deletitions
find "$DIR" -name "*.md" -exec sed -i 's/\\tabularnewline//g' {} +
find "$DIR" -name "*.md" -exec sed -i 's/\\hline//g' {} +

echo "Replacements completed."