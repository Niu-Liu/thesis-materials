fin = open('AMS262.dat','r')
lin = fin.readlines()
fou = open('AMS262.cat','w')
for i in range(len(lin)):
	print >> fou, lin[i][22:30]
fou.close()
fin.close()
