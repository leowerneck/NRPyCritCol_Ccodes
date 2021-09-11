# First compile the code
echo "Compiling NRPyCritCol..."
make
echo "Done!"

# Set common initial data parameters
pulse_center=0
pulse_width=1
Nr=200001
R_max=70
outputfile=outputSFC.txt

# Generate weak initial data
echo "Generating weak initial data..."
pulse_amplitude=0.30332394090
python ScalarField_Initial_Data.py $outputfile $pulse_amplitude $pulse_center $pulse_width $Nr $R_max
cp $outputfile ID_weak.txt
echo "Done! Saved a copy of initial data file to ID_weak.txt."

# Run the code
echo "Running NRPyCritCol for weak initial data..."
./NRPyCritCol 320 2 2 64 0.2 2.293577981651376 0.0

# Move output file
mv out_central_values.dat out_weak.dat
echo "Done! Central values stored in out_weak.dat."

# Generate strong initial data
echo "Generating strong initial data..."
pulse_amplitude=0.30332394095
python ScalarField_Initial_Data.py $outputfile $pulse_amplitude $pulse_center $pulse_width $Nr $R_max
cp $outputfile ID_strong.txt
echo "Done! Saved a copy of initial data file to ID_strong.txt."

# Run the code
echo "Running NRPyCritCol for strong initial data..."
./NRPyCritCol 320 2 2 64 0.2 2.293577981651376 0.0

# Move output file
mv out_central_values.dat out_strong.dat
echo "Done! Central values stored in out_strong.dat."

# Generate the plot
echo "Generating lapse self-similar behaviour plot..."
python generate_plot.py
echo "Done! Saved plot to file lapse_self_similarity.png."

# All done!
