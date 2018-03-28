import os

datadir = '/mnt/75/16s'#/media/vd1/171213_16s'
file_name_list = []
file_name_1 = ''
file_name_2 = ''
dirname = ''

def get_file_list(dir):
	for f in os.listdir(dir):
		if 'XK' in f:
			file_name_list.append(f)
	#print(file_name_list) 
	return file_name_list


def get_operation_file_names(file):
	if '_R1' in file:
		file_name_1 = file
		file_name_2 = file.split('_R1')[0] + '_R2' + file.split('_R')[1]
		dirname = file.split('_R')[0]
		print file_name_1, file_name_2,dirname
		return file_name_1, file_name_2, dirname
	

file_name_list = get_file_list(datadir)

for file in file_name_list:
	file_name_1, file_name_2, dirname = get_operation_file_names(file)
	cmd = 'bash ./16S-processing.sh %s %s %s' % (file_name_1,  file_name_2, dirname)
	os.makedirs(dirname)
	print cmd
	os.system(cmd)
	break


