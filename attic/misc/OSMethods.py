#!/usr/bin/env/python
##============================================================

import os, glob
import shutil

##============================================================

## Copy files with hard-coded string to public folder

def main():
	
	for string in ["CFHTLS_*_sci.png"]:	
		for image in file_seek("../Data/",string):
			shutil.copy(image,"/home/users2/cs3006/public_html/PSF_Eyeball")	
	return

##============================================================

## Walks through directories to find file of choice
	## ext = "*.fits" for example
def file_seek(top_dir, ext):
	
	file_list = []
	
	## Walk through the directories, starting with top_dir
	for root,dirs,files in os.walk(top_dir):
		## Find .ext file
		for my_file in glob.iglob( os.path.join( root, ext )):
			## Don't use homogenisation kernels
			if my_file[-9:-5]!="homo":	file_list += [my_file]
	return file_list
	
##============================================================
## Delete all files with given (hard-coded) strings
	## From Python line
def del_string():	
	for string in ["*_.png"]:#"*_cBB.png","*_cHE.png","*_br8.png"
		for image in file_seek("../Data/",string):
			os.remove(image)
			print "OSScript.py: del_string: deleted file",image
	return

##============================================================
if __name__=="__main__":
	main()
