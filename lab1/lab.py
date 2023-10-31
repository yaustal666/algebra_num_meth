from random import randint

with open("matrix.txt", "w", encoding="utf-8") as f:

    for i in range(2):
        for i in range(9):
            a = randint(1, 10)
            f.write(str(a))
            f.write(" ")
        f.write("\n")

# with open("matrix.txt", "w", encoding="utf-8") as f:

#     for i in range(2):
#         for i in range(625):
#             a = randint(1, 10)
#             f.write(str(a))
#             f.write(" ")
#         f.write("\n")