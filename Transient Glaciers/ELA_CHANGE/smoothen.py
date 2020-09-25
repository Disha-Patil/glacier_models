with open("bed.txt") as f: b_array = [[int(j[0]),int(j[1]),float(j[2])] for j in [k.split(" ") for k in f.readlines()]]
xmax,ymax = b_array[-1][0]+1,b_array[0][1]+1

bed = [[0 for j in range(ymax)] for i in range(xmax)]
smoothbed = bed
mask = [[0 for j in range(ymax)] for i in range(xmax)]

for t in b_array:
    [i,j,k] = t
    bed[i][j] = k
    if(k): mask[i][j] = 1

for n in range(1):
    blank = [[0 for j in range(ymax)] for i in range(xmax)]
    iv = [0, 1, -1, 0, 0, 1, 1, -1, -1]
    jv = [0, 0, 0, 1, -1, 1, -1, 1, -1]
    w  = [4,2,2,2,2,1,1,1,1]
    for i in range(xmax):
        for j in range(ymax):
            if(mask[i][j]):
                for k in range(9): blank[i][j] += smoothbed[i+iv[k]][j+jv[k]]*w[k]
                total = 0
                for k in range(9): total += w[k]*mask[i+iv[k]][j+jv[k]]
                blank[i][j] /= total
    smoothbed = [row[:] for row in blank]

lm = open( "smoothbed.txt", 'w')
for j in range(ymax-1,-1,-1):
    for i in range(xmax):
        lm.write('%d %d %f\n' %(i,j,smoothbed[i][j]))
lm.close()
