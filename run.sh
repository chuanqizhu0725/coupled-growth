g++ -I /opt/X11/include -L /opt/X11/lib -lX11 \
*.cpp -o main && rm -f data/phi/*.vtk data/phi/*.csv figures/phi/*.png data/temp/*.csv figures/temp/*.png && ./main && python plot1d.py && rm main