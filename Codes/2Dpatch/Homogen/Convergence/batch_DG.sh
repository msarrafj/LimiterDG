#! /bin/bash
mesh="2 4 8 16 32"
echo -ne "nx \t (dofs) \t L2_Sw \t L2_Pw\t H1_Sw\t H1_Pw\n" > log_DG.dat
for i in $mesh
    do
    echo "Solving for nx = $i"
    python DG.py ${i} > tmp_DG.dat
    paste <(echo -ne "${i}") <(awk '$1=="dofs" {printf "%d\n" ,$2}' tmp_DG.dat|tail -n 1) <(awk '$1=="L2_error_of_s" {printf "%e\n" ,$2}' tmp_DG.dat|tail -n 1) <(awk '$1=="L2_error_of_p" {printf "%e\n" ,$2}' tmp_DG.dat|tail -n 1) <(awk '$1=="H1_error_of_s" {printf "%e\n" ,$2}' tmp_DG.dat|tail -n 1) <(awk '$1=="H1_error_of_p" {printf "%e\n" ,$2}' tmp_DG.dat|tail -n 1) >>log_DG.dat
done
rm tmp_DG.dat
