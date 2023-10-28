from random import randint

with open("matrix.txt", "w", encoding="utf-8") as f:

    for i in range(2):
        for i in range(16):
            a = randint(1, 6)
            f.write(str(a))
            f.write(" ")
        f.write("\n")