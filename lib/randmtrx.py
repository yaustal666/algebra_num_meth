from random import randint
import pandas as pd
msize = 25

mtrx = []
nums = []
idline = [1] * msize

for i in range(msize):
    nums.append(randint(0, 999) / 1000)

# random Vandermond matrix
with open("matrix.txt", "w", encoding="utf-8") as f:
    
    k = 0
    for i in range (msize):
        tmp = []
        for j in nums:
            a = pow(j, k)
            f.write(str(a))
            f.write(" ")
            tmp.append(a)
        
        mtrx.append(tmp)
        k += 1
        f.write("\n")
matrix = pd.DataFrame(mtrx)
matrix.to_csv("mtrx_csv.csv")


# random matrix
# mtrx = []
# with open("matrix.txt", "w", encoding="utf-8") as f:

#     for i in range(msize):
#         tmp = []
#         for i in range(msize):
#             a = randint(1, 999) / 1000
#             tmp.append(a)
#             f.write(str(a))
#             f.write(" ")
#         mtrx.append(tmp)
#         f.write("\n")

# matrix = pd.DataFrame(mtrx)
# matrix.to_csv("mtrx_csv.csv")