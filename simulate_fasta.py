import sys
import random

n = int(sys.argv[1])
m = int(sys.argv[2])
d = float(sys.argv[3])

num_base_map = {0:"A", 1:"T", 2:"G", 3:"C"}
mutation_map = {"A":"TGC", "T":"AGC", "G":"ATC", "C":"ATG"}
S = ""
for i in range(m):
    S += num_base_map[random.randint(0, 3)]

print("> seq" + str(0))
print(S)
print("")

for l in range(1, n):
    result = ""
    for i in range(m):
        if random.randint(1, 100) <= int(d*100):
            result += mutation_map[S[i]][random.randint(0, 2)]
        else:
            result += S[i]
    print("> seq" + str(l))
    print(result)
    print("")