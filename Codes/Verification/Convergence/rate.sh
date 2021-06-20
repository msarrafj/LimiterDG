#! /bin/bash
awk 'BEGIN{OFS="\t"}NR==1{print $0,"rateL2_Sw"}NR>1{print $0,-(log($3)-b)/(log($1)-a),-(log($4)-c)/(log($1)-a),-(log($5)-d)/(log($1)-a),-(log($6)-e)/(log($1)-a)}{a=log($1);b=log($3);c=log($4);d=log($5);e=log($6)}' $1

