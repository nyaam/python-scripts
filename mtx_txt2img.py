#!/usr/bin/env python

"""
Convert protein sequence to image
bulid time: 25MAR2018
author: N. J. Cheung
email: nyaam.ch@gmail.com
license: MIT
"""
import sys
import getopt
import numpy
import math
from PIL import ImageFont, Image, ImageDraw


"""
>3GTI:A|PDBID|CHAIN|SEQUENCE
"""
# seq     = 'HHHHHHGLVPRGSHMTAQTVTGAVAAAQLGATLPHEHVIFGYPGYAGDVTLGPFDHAAALASCTETARALLARGIQTVVD \
#            ATPNGCGRNPAFLREVSEATGLQILCATGFYYEGGGATTYFKFRASLGDAESEIYEMMRTEVTEGIAGTGIRAGVIKLAS \
#            SRDAITPYEQLFFRAAARVQRETGVPIITHTQEGQQGPQQAELLTSLGADPARIMIGHMDGNTDPAYHRETLRHGVSIAF \
#            DRIGLQGLVGTPTDAERLSVLTTLLGEGYADRLLLSHDSIWHWLGRPPAIPEAALPAVKDWHPLHISDDILPDLRRRGIT \
#            EEQVGQMTVGNPARLFG'
# seq  =  seq.replace(" ", "")
# YYYYYYYYYYY
# AAAAAYAAAAA
# YYYYYYYYYYY

# AAAAAYAAAAA
# YYYYYYYYYYY
# AAAAAYAAAAA

name    = 'text2img'
version = 1.0
output_filename = 'default.png'

AA = 20
wsz = 5

w = 28
h = w
# TODO:
# Better way to convert AA to data matrix?????
def usage():
	print("NAME")
	print("\t%s" % name)
	print("VERSION")
	print("\t%.1f" % version)
	print("SYNOPSIS")
	print("\tpython3 %s.py [-h] [-t text -m matrix]" % (name))
	print("DESCRIPTION")
	print("\t%s.py converts both text and data matrix to an image" % name)
	print("OPTIONS")
	print("\t-h --help \t Display help info.")
	print("\t-t --text \t Text in ASCII")
	print("\t-m --matrix \t Data matrix")
	print("\t-o --output \t Output image in *.png")
	print("\t-version \t Version")
	print("AUTHOR")
	print("\tN. J. Cheung <nyaam.ch@gmail.com>")
	sys.exit(0)

# Slide window over the given text
def txt_slide(txt_fid):
	name_ = txt_fid.split('.')[0]
	seq = read_txt(txt_fid)
	l = len(seq)
	for i in range(l):
		s = ''
		if i < wsz:
			for j in range(wsz):
				s += seq[i]
			# s += '-'+seq[i]+'-'#seq[i]
			s += seq[i]#seq[i]
			s += seq[i+1:i+wsz+1]
			# print(s)
		if i > wsz-1 and i < l-wsz:
			s = seq[i-wsz:i]
			# s += '-'+seq[i]+'-'
			s += seq[i]
			s += seq[i+1:i+wsz+1]
			# print(s)
		if i > l-wsz-1:
			# s = ''
			s += seq[i-wsz:i]
			# s += '-'+seq[i]+'-'
			s += seq[i] 
			s += seq[i+1:]
			s += seq[i]*(i+wsz-l+1)
		mtx_seq = []#numpy.array([])
		for j in s:
			mtx_seq.append((ord(j)-65)/AA)
		mtx_seq = numpy.asarray(mtx_seq)
		mtx_seq_T = mtx_seq.reshape((mtx_seq.shape[0],1))

		width  = 23
		height = wsz*2+1
		for j in range(w-width):
			if j < 1:
				mtx = mtx_seq_T #numpy.concatenate((mtx_T, mtx), axis = 1)
			else:
				mtx = numpy.concatenate((mtx, mtx_seq_T), axis = 1)

		# Normalization of AA
		# row_sums = mtx.sum(axis = 1)
		# mtx = mtx / row_sums[:, numpy.newaxis]

		# Features: PSSM+SS+SA
		dat = numpy.random.random((height, width))
		# Normalization of features
		# row_sums = dat.sum(axis = 1)
		# dat = dat / row_sums[:, numpy.newaxis]

		mtx = numpy.concatenate((mtx, dat), axis = 1)
		# Normalization of all
		row_sums = mtx.sum(axis = 1)
		mtx = mtx / row_sums[:, numpy.newaxis]
		text2img(name_, s, mtx)
	# text2img(name_, 'ABCDEFGHIJK', mtx)
	# text2img(name_, 'HHHHHHHHHHH', mtx)
	# text2img(name_, 'MMMMMMMMMMM', mtx)
	return s

def text2img(name_, txt, mtx):
	l = len(txt)
	idx = math.floor(l/2)  # Index of current amino acid
	aa = txt[idx]
	one = numpy.ones((w-mtx.shape[0], 1))*mtx[idx]
	mtx_ = numpy.concatenate((mtx[0:idx].T, one.T), axis = 1)
	mtx_ = numpy.concatenate((mtx_, mtx[idx:].T), axis = 1)
	mtx = mtx_.T
	
	name_img = "img/" + name_ + '-' + aa + ".png"
	# Rescale to 0-255 and convert to uint8
	dat = (255.0 / mtx.max() * (mtx - mtx.min())).astype(numpy.uint8)
	im = Image.fromarray(dat)
	# dat = numpy.transpose(dat, (1, 0))
	# im.show()
	im.save(name_img)

# Read text from file
def read_txt(input_file):
	txt = ''
	with open(input_file, "r") as fid:
		for line in fid:
			line = line[:-1]
			if line[0] != '>':
				txt += line
	return txt

# Read matrix
def read_mtx(mtx_fid):
	# txt = ''
	with open(mtx_fid, "r") as fid:
		for line in fid:
			line = line[:-1]
	# return txt
# Write matrix
# def write_mtx(img_fid, mtx):
# 	with open (img_fid, "w") as fid:
# 		for line in mtx:
# 			# print(line)
# 			fid.write(line + "\n")

def main():
	if(not len(sys.argv[1:])):
		usage()
	try:
		opts, args = getopt.getopt(sys.argv[1:], "t:m:o:hv", ["text=", "matrix=", "output=", "help", "version"])
	except getopt.GetoptError as error:
		print(str(error))
		usage()
	if len(opts) < 2:
		usage()
	for opt, arg in opts:
		if   opt in ("-h", "--help"):
			usage()
		elif opt in ("-t", "--text"):
			txt_fid = arg
		elif opt in ("-m", "--matrix"):
			mtx_fid = arg
		elif opt in ("-o", "--output"):
			img_fid = arg
		elif opt == '--version':
			print("%s version %.1f" % (name, version))
			sys.exit(0)
		else:
			assert False, "Unrecognized option"
			sys.exit(-1)

	txt_slide(txt_fid)

if __name__ ==  "__main__()":
	main()

main()


