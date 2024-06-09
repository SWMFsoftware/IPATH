
#!/bin/bash



mpif90 -O3 parallel_wrapper.f transport_back_1.09.f -o proton_80.out
#mpif90 -O3 parallel_wrapper.f iron_1.09.f -o iron_80.out -lhwloc
#mpif90 -O3 parallel_wrapper.f cno_1.09.f -o cno_80.out -lhwloc

#mkdir $root_dir/proton_80_1.09
#mkdir $root_dir/iron_80_1.09
#mkdir $root_dir/cno_80_1.09

gfortran combine.f -o combine.out

mv proton_80.out proton_80_1.09
#mv iron_80.out $root_dir/iron_80_1.09
#mv cno_80.out $root_dir/cno_80_1.09

cp combine.out proton_80_1.09
#cp combine60.out scpt_iron $root_dir/iron_80_1.09
#cp combine60.out scpt_cno $root_dir/cno_80_1.09

cp ../CME/path_output/dist_at_shock.dat proton_80_1.09
cp ../CME/path_output/shock_momenta.dat proton_80_1.09
cp ../CME/path_output/esc_distr_dn.dat proton_80_1.09
cp ../CME/path_output/all_shell_bndy.dat proton_80_1.09
cp ../CME/path_output/dist_all_shl.dat proton_80_1.09

cd ./proton_80_1.09
rm fp*
rm RawData_*
mpirun -np 4 ./proton_80.out
./combine.out
rm Raw*
