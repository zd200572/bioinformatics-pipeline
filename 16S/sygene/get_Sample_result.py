import os

os.makedir('p_')
os.makedir('g_')

for dir_name in os.listdir('.'):
	if 'XK' in dir_name:
		cmd1 = 'cp dir_name/result/sum_taxa/otu_table3_L2.txt p_/%sotu_table3_L2.txt' % dir_name
		cmd2 = 'cp dir_name/result/sum_taxa/otu_table3_L6.txt g_/%sotu_table3_L6.txt' % dir_name
		os.system(cmd1)
		os.system(cmd2)
