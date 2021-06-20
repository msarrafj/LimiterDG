#! /bin/bash
mesh="2 4 8 16 32"
echo -ne "nx \t (dofs) \t L2_Sw \t L2_pw\t H1_Sw\t H1_pw\n" > log_DG_avg.dat
for i in $mesh
    do
    echo "Solving for nx = $i"
    python DG_avg.py ${i} > tmp_avg.dat
    paste <(echo -ne "${i}") <(awk '$1=="dofs" {printf "%d\n" ,$2}' tmp_avg.dat|tail -n 1) <(awk '$1=="L2_error_of_s" {printf "%e\n" ,$2}' tmp_avg.dat|tail -n 1) <(awk '$1=="L2_error_of_p" {printf "%e\n" ,$2}' tmp_avg.dat|tail -n 1) <(awk '$1=="H1_error_of_s" {printf "%e\n" ,$2}' tmp_avg.dat|tail -n 1) <(awk '$1=="H1_error_of_p" {printf "%e\n" ,$2}' tmp_avg.dat|tail -n 1) >>log_DG_avg.dat
done
rm tmp_avg.dat
