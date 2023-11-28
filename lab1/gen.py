from random import randint

mtrx = []
mwidth = 3
mheight = 3

with open("matrix.txt", "w", encoding="utf-8") as f:
    f.write(str(mwidth))
    f.write(" ")
    f.write(str(mheight))
    f.write("\n")
    for i in range(mheight):


        tmp = []
        for i in range(mwidth):
            a = randint(1, 20)
            tmp.append(a)
            f.write(str(a))
            f.write(" ")
        mtrx.append(tmp)
        f.write("\n")