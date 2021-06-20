SIZE = 5;
Tol = 0.1;
p1x = 20;
p1y = 35;
p2x = 75;
p2y = 40;
p3x = 25;
p3y = 15;
p5x = 60;
p5y = 65;
p4x = 60;
p4y = 50;
p8x = 42;
p8y = 80;
p6x = 40;
p6y = 55;
p7x = 65;
p7y = 85;
Point(1) = {0, 0, 0, SIZE};
Point(2) = {100, 0, 0, SIZE};
Point(3) = {100, 100, 0, SIZE};
Point(4) = {0, 100, 0, SIZE};
// Cracks 1: point 1 to point 2
Point(5) = {p1x, p1y, 0, SIZE};
Point(6) = {p2x, p2y, 0, SIZE};
/* Point(7) = {0.5*(p1x+p2x)+Tol, 0.5*(p1y+p2y)-Tol, 0, SIZE}; */
/* Point(8) = {0.5*(p1x+p2x)-Tol, 0.5*(p1y+p2y)+Tol, 0, SIZE}; */
//+
// Cracks 2: point 3 to point 5
Point(9) = {p3x, p3y, 0, SIZE};
Point(10) = {p5x, p5y, 0, SIZE};
/* Point(11) = {0.5*(p3x+p5x)+Tol, 0.5*(p3y+p5y)-Tol, 0, SIZE}; */
/* Point(12) = {0.5*(p3x+p5x)-Tol, 0.5*(p3y+p5y)+Tol, 0, SIZE}; */
//+
// Cracks 3: point 4 to point 8
Point(13) = {p4x, p4y, 0, SIZE};
Point(14) = {p8x, p8y, 0, SIZE};
/* Point(15) = {0.5*(p4x+p8x)+Tol, 0.5*(p4y+p8y)+Tol, 0, SIZE}; */
/* Point(16) = {0.5*(p4x+p8x)-Tol, 0.5*(p4y+p8y)-Tol, 0, SIZE}; */
//+
// Cracks 4: point 6 to point 7
Point(17) = {p6x, p6y, 0, SIZE};
Point(18) = {p7x, p7y, 0, SIZE};
/* Point(19) = {0.5*(p6x+p7x)+Tol, 0.5*(p6y+p7y)-Tol, 0, SIZE}; */
/* Point(20) = {0.5*(p6x+p7x)-Tol, 0.5*(p6y+p7y)+Tol, 0, SIZE}; */
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Point(21) = {40.081, 36.736, 0, 0.075};
//+
Point(22) = {40.397, 36.759, 0, 0.075};
//+
Point(23) = {40.523, 36.947, 0, 0.075};
//+
Point(24) = {40.184, 36.933, 0, 0.075};
//+
Point(25) = {55.129, 57.973, 0, 0.075};
//+
Point(26) = {55.225, 58.109, 0, 0.075};
//+
Point(27) = {55.177, 58.187, 0, 0.075};
//+
Point(28) = {55.082, 58.049, 0, 0.075};
//+
Point(29) = {50.021, 66.85, 0, 0.075};
//+
Point(30) = {49.907, 67.064, 0, 0.075};
//+
Point(31) = {49.733, 66.856, 0, 0.075};
//+
Point(32) = {49.867, 66.656, 0, 0.075};
//+
Line(5) = {5, 21};
//+
Line(6) = {21, 9};
//+
Line(7) = {9, 22};
//+
Line(8) = {22, 6};
//+
Line(9) = {6, 23};
//+
Line(10) = {23, 25};
//+
Line(11) = {25, 13};
//+
Line(12) = {13, 26};
//+
Line(13) = {26, 10};
//+
Line(14) = {10, 27};
//+
Line(15) = {27, 29};
//+
Line(16) = {29, 18};
//+
Line(17) = {18, 30};
//+
Line(18) = {30, 14};
//+
Line(19) = {14, 31};
//+
Line(20) = {31, 17};
//+
Line(21) = {17, 32};
//+
Line(22) = {32, 28};
//+
Line(23) = {28, 24};
//+
Line(24) = {24, 5};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {18, 19, 20, 21, 22, 23, 24, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve(1) = {4};
//+
Physical Curve(2) = {2};
//+
Physical Surface(1) = {1};
//+
Transfinite Curve {9, 8} = 17 Using Progression 1;
//+
Transfinite Curve {17, 16, 10, 10, 23, 6, 7} = 12 Using Progression 1;
//+
Transfinite Curve {18, 18, 19, 20, 20, 21, 5, 5, 24} = 9 Using Progression 1;
//+
Transfinite Curve {15, 15, 22, 14, 13, 12, 11} = 7 Using Progression 1;
//+
Transfinite Curve {4, 3, 2, 1} = 13 Using Progression 1;
