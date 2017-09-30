import os
for  i in range(205, 210):
	print(i)
	#os.system("ascp -QT")
	cmd = ' ascp -QT -l 6M -k1 -i asperaweb_id_dsa.openssh anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP017/SRP017311/SRR620' + '%s' % i +'/SRR620' + '%s.sra' % i + ' d:\\'
	print(cmd)
	os.system(cmd)
	#break