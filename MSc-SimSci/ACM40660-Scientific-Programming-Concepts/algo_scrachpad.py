import numpy

//extract elements

a = np.matrix('1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16')

b = np.matrix('0 0 0; 0 0 0; 0 0 0')

r = 0
c = 0
b_r = 0
b_c = 0
for i in [x for x in range(0, 4) if x != r]:
    b_c = 0
    for j in [x for x in range(0, 4) if x != c]:
        b[b_r, b_c] = a[i,j]
        b_c = b_c + 1
    b_r = b_r  + 1
        
b