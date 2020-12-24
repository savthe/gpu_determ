f1 = open("true_results.out", 'r')
f2 = open("results.out", 'r')
lines1 = f1.readlines()
lines2 = f2.readlines()
assert len(lines1) == len(lines2)

result = [0., 0., 0., 0., 0.]
for i in range(0, len(lines1)):
	x = lines1[i].split(' ')
	y = lines2[i].split(' ')
	for k in range(0, len(result)):
		result[k] = result[k] + (float(x[4+k]) - float(y[4+k]))**2

for k in range(0, len(result)):
	print((result[k] / len(lines1))**(.5), end = '\t')
print()

