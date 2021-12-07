res = 50
for j in range(res):
	if j <= res / 2:
		J = 10. + float(j) / float(res/2) * 90.
	else:
		J = float(j-res/2) / float(res/2-1) * 900. + 100.
	print(j,J)
