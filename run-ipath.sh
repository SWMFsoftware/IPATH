#!/bin/sh

# Generate a steady background solar wind
cd ./backsw
rm zr*
rm xdzeus36
csh -v dzeus36.s
./xdzeus36

# Generate a cme at given location and time
cd ../CME
rm zhto*
rm zr*
rm xdzeus36
cp ../backsw/zr004TW .
csh -v dzeus36.s
./xdzeus36 <input

# Transport
cd ../src
chmod u+x run_trspt.sh
./run_trspt.sh

# Plot CME figures and generate a movie
cd ../plotting
python swprofile.py
python plot-density_simple.py
python shock_param_overview.py
python time_intens_new.py
python ESP_new.py
convert -delay 10 CME-Density*.png -loop 0 CME-Density_animated.gif
mv CME-Density_animated.gif CME-Density_-`date +%Y%m%d%H%M`.gif
mkdir `date +%Y%m%d%H%M`
cp ./*.png ./`date +%Y%m%d%H%M`/.
cp ./*.gif ./`date +%Y%m%d%H%M`/.
rm *.png
rm *.gif

