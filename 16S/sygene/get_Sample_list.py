import os

fout = open('Samples_list.csv','w')
fout.write("sample-id,absolute-filepath,direction" + '\n')
path = '/home/zjd/16s/test'
for i in os.listdir(path):
	if 'R1.fastq.gz' in i:
		sample_id = i.split('__R1')[0]
		fout.write(str(sample_id) + ',' + path + '/' + sample_id + '_1.fastq.gz' +',' + 'forward' + '\n')
		fout.write(str(sample_id) + ',' + path + '/' + sample_id + '_2.fastq.gz' +',' + 'reverse' + '\n')
fout.close()
