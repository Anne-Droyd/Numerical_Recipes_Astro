#!/bin/bash

echo "Run handin kirwan.tex)!"

echo "Clearing/creating the plotting directory"
if [ ! -d "plots" ]; then
  mkdir plots
fi
rm -rf plots/*


# Script that returns a plot
echo "Running the hand in for Kirwan: "
python3 hand_ins/Hand_in_1/code/Hand_in_1_code.py

echo "Generating the pdf"
cd hand_ins/Hand_in_1/code
pdflatex Hand_in_1_pdf.tex