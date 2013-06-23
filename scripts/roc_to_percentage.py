import sys

lines = sys.stdin.readlines()

total_c0 = 0
total_c1 = 1
last_cc0 = -1

for line in lines:
    if line == None:
        break;

    tokens = line.split()
    #print tokens

    c1 = int(tokens[1])
    c0 = int(tokens[3])
    cc0 = int(tokens[4])

    total_c0 += c0
    total_c1 += c1
    last_cc0 = cc0


print "%d/%d(%d) = %0.2lf\n" % (total_c0-total_c1, total_c0, last_cc0, 100.0*(total_c0-total_c1)/total_c0)
