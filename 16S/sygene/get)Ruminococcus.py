with open('genus.normalization.profile') as f:
	for line in f:
		if line.strip().split('	')[0] == 'Ruminococcus':
			for l in line.strip().split('	'):
				if l > str(0.2):
					print(l)
