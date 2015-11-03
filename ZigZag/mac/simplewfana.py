# To run script type python path/simplewfana.py nnbar-rawdigit-file numi-truth-file numi-rawdigit-file

import sys,os

if len(sys.argv) < 2:
	msg  = '\n'
	msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
	msg += '\n'
	sys.stderr.write(msg)
	sys.exit(1)

from ROOT import gSystem
#gSystem.Load("libSimpleWFAna")
from ROOT import larlite as fmwk

fname = "ttree_output.root"

T = 10

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# add input files - nnbar rawdigit file
my_proc.add_input_file(sys.argv[1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("SimpleWFAna_output.root");

ana_unit = fmwk.SimpleWFAna(T,fname)

my_proc.add_process(ana_unit)

print
print  "Finished configuring ana_processor. Start event loop!"
print
		
# Let's run it.
my_proc.run()

# done!
print 
print "Finished running ana_processor event loop!"
print

sys.exit(0)

