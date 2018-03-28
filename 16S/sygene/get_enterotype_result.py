import os

dic = {}

for fi in os.listdir('p_'):
	#print(fi)
	di = {}
	if 'XK' not in fi:
		continue
	with open('p_/%s' % fi) as f:
		for line in f:
			if 'Firmicutes' in line:
				F = line.strip().split('	')[0]
				di[F] = line.strip().split('	')[1]
			elif 'Bacteroidetes' in line:
				B = line.strip().split('	')[0]
				di[B] = line.strip().split('	')[1]
		if di[F] > di[B]:
			print('F:', di[F], 'B:', di[B], 'ETF')
	fi2 = str(fi.split('_')[0] + '_table3_L6.txt')
	print(fi.split('otu')[0])
	if fi2 in os.listdir('g_'):
		with open('g_/%s' % fi2) as f2:
			for line in f2:
				if 'Bacteroides' in line:
					Bacteroides = 'Bacteroides'
					di[Bacteroides] = line.strip().split('	')[1]
				if "Prevotella" in line:
					Prevotella = 'Prevotella'
					di[Prevotella] = line.strip().split('	')[1]
				if 'Ruminococcus' in line:
					Ruminococcus = 'Ruminococcus'
					di[Ruminococcus] = line.strip().split('	')[1]
			if 'Prevotella' not in di .keys():
				print('Bacteroides:', di[Bacteroides], 'Ruminococcus:', di[Ruminococcus], 'ETB')
				continue
			if di[Bacteroides] > di[Prevotella]:
				print('Bacteroides:', di[Bacteroides], 'Prevotella:', di[Prevotella], 'Ruminococcus:', di[Ruminococcus], 'ETB')
			else:
				print('Bacteroides:', di[Bacteroides], 'Prevotella:', di[Prevotella], 'Ruminococcus:', di[Ruminococcus], 'ETP')
	dic[fi.split('_')[0]] = di
	#break