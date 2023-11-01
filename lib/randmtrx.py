from random import randint

msize = 25

with open("matrix.txt", "w", encoding="utf-8") as f:

    for i in range(msize):
        for i in range(msize):
            a = randint(1, 999) / 1000
            f.write(str(a))
            f.write(" ")
        f.write("\n")
