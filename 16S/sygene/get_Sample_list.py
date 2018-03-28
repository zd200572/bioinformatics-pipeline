import os

fout = open('Samples_list.csv','w')
fout.write("sample-id,absolute-filepath,direction" + '\n')
for i in os.listdir('.'):
	if 'R1.fq.gz' in i:
		sample_id = i.split(_R1)[0]
		fout.write(str(sample_id) + ',' + str(os.path.getcwd()) + '/' + sample-id + '_1.fq.gz' +',' + 'forward' + '\n')
		fout.write(str(sample_id) + ',' + str(os.path.getcwd()) + '/' + sample-id + '_2.fq.gz' +',' + 'reverse' + '\n')
fout.close()