import os, shutil

with open('./init_cond.txt', 'r') as fin:
	mean = int( fin.readline() )

for name in os.listdir():
	if name.split( sep = '.' )[-1] == 'out':
		new_name = './old/' + name.split( sep = '.' )[0] + '_{0}'.format(mean+1) + '.out'
		shutil.move(name,new_name)
