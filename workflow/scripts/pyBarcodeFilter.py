#!/usr/bin/python
# updated by vincent rocher to fit snakemake and change the spacer of BarcodeSplitter (by "_")
__authors__		= ["Sander Granneman"]
__copyright__	= "Copyright 2019"
__version__		= "3.1"
__credits__		= ["Sander Granneman","Grzegorz Kudla"]
__maintainer__	= "Sander Granneman"
__email__		= "sgrannem@ed.ac.uk"
__status__		= "Production"

##################################################################################
#
#	pyBarcodeFilter.py
#
#
#	Copyright (c) Sander Granneman and Grzegorz Kudla 2019
#	
#	Permission is hereby granted, free of charge, to any person obtaining a copy
#	of this software and associated documentation files (the "Software"), to deal
#	in the Software without restriction, including without limitation the rights
#	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#	copies of the Software, and to permit persons to whom the Software is
#	furnished to do so, subject to the following conditions:
#	
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#	
#	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#	THE SOFTWARE.
#
##################################################################################

import sys
import re
import time
from optparse import *
from pyCRAC.Classes.Barcodes import *
import shutil

data = BarcodeSplitter("fastq",snakemake.input["fastq"],"",snakemake.input["barcode"],allowedmismatches=0,gzip=False,keepbarcode=False,search="forward")
data.spacer="_"
data.processFastQFile()
data.printRandBarcodeStats()
data.printBarcodeStats()
