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
"""
print('Enter a cut (1=hit num, 2=TDC std, 3=ADC amp, 4=integrated WF, 5=TDC iqr):')
while True:
	try:
		option = int(raw_input());
		break
	except(option not in [1,2,3,4,5]):
		print('Must be an integer between 1 and 5')

print('Enter a plane(s) to plot energy dist (1=U, 2=V, 3=Y, 4=UV, 5=UY, 6=VY, 7=UVY):')
while True:
        try:
                plane = int(raw_input());
                break
        except(option not in [1,2,3,4,5,6,7]):
                print('Must be an integer between 1 and 7')
"""
plane = 4
for y in [0,1,4,2,3]:
	if y == 0:
		option = 1
		op_name = "HitNum"
	if y == 1:
		option = 2
		op_name = "TDCstd"
	if y == 2:
		option = 3
		op_name = "ADCamp"
	if y == 3:
		option = 4
		op_name = "WFint"
	if y == 4:
		option = 5
		op_name = "TDCiqr"
	# Loop over different tolerances changing file names each time
	for x in range (3,4):
		if x == 0:
			T = 4
			fname = op_name+"-nnbar4.root"
			fname1 = op_name+"-numi4.root"
			tname = op_name+"-nnbar4.txt"
			tname1 = op_name+"-numi4.txt"
			print
			print(op_name+", T = "+str(T))
			print
		if x == 1:
			T = 6
			fname = op_name+"-nnbar6.root"
			fname1 = op_name+"-numi6.root"
			tname = op_name+"-nnbar6.txt"
			tname1 = op_name+"-numi6.txt"
                        print
                        print(op_name+", T = "+str(T))
                        print
		if x == 2:
			T = 8
			fname = op_name+"-nnbar8.root"
			fname1 = op_name+"-numi8.root"
			tname = op_name+"-nnbar8.txt"
			tname1 = op_name+"-numi8.txt"
                        print
                        print(op_name+", T = "+str(T))
                        print
		if x == 3:
			T = 10
			fname = op_name+"-nnbar10.root"
			fname1 = op_name+"-numi10.root"
			tname = op_name+"-nnbar10.txt"
			tname1 = op_name+"-numi10.txt"
                        print
                        print(op_name+", T = "+str(T))
                        print
		if x == 4:
			T = 12
			fname = op_name+"-nnbar12.root"
			fname1 = op_name+"-numi12.root"
			tname = op_name+"-nnbar12.txt"
			tname1 = op_name+"-numi12.txt"
                        print
                        print(op_name+", T = "+str(T))
                        print
		if x == 5:
			T = 14
			fname = op_name+"-nnbar14.root"
			fname1 = op_name+"-numi14.root"
			tname = op_name+"-nnbar14.txt"
			tname1 = op_name+"-numi14.txt"
                        print
                        print(op_name+", T = "+str(T))
                        print

		# Create ana_processor instance
		my_proc = fmwk.ana_processor()
	
		# add input files - nnbar rawdigit file
		my_proc.add_input_file(sys.argv[1])
	
		# Specify IO mode
		my_proc.set_io_mode(fmwk.storage_manager.kREAD)
	
		ana_unit = fmwk.SimpleWFAna(T,fname,tname,option,plane)
	
		my_proc.add_process(ana_unit)
	
		print
		print  "Finished configuring ana_processor. Start event loop!"
		print
		
		# Let's run it.
		my_proc.run()
		# Get cut limits from nnbar file
		t_min = ana_unit.GetTmin()
		t_max = ana_unit.GetTmax()
		u_min = ana_unit.GetUmin()
		u_max = ana_unit.GetUmax()
		v_min = ana_unit.GetVmin()
		v_max = ana_unit.GetVmax()
		y_min = ana_unit.GetYmin()
		y_max = ana_unit.GetYmax()
		# done!
		print 
		print "Finished running ana_processor event loop!"
		print

		# Create ana_processor instance
		my_proc1 = fmwk.ana_processor()

		# add input files - numi truth file
		my_proc1.add_input_file(sys.argv[2])

		# Specify IO mode
		my_proc1.set_io_mode(fmwk.storage_manager.kREAD)

		ana_unit1 = fmwk.SimpleWFAna(T,"delete.root","delete.txt",option,plane)

		my_proc1.add_process(ana_unit1)

		print
		print  "Finished configuring ana_processor. Start event loop!"
		print

		# Let's run it.
		my_proc1.run()
		# Get interaction types
		typ = ana_unit1.GetType()
		# Get neutrino energies
		qsq = ana_unit1.GetQsq()
		# done!
		print
		print "Finished running ana_processor event loop!"
		print

		# Create ana_processor instance
		my_proc2 = fmwk.ana_processor()

		# add input files
		my_proc2.add_input_file(sys.argv[3])
	
		# Specify IO mode
		my_proc2.set_io_mode(fmwk.storage_manager.kREAD)
	
		ana_unit2 = fmwk.SimpleWFAna(T,fname1,tname1,option,plane)
	
		# Set cut limits, interaction types and neutrino energies
		ana_unit2.SetTmin(t_min)
		ana_unit2.SetTmax(t_max)
		ana_unit2.SetUmin(u_min)
		ana_unit2.SetUmax(u_max)
		ana_unit2.SetVmin(v_min)
		ana_unit2.SetVmax(v_max)
		ana_unit2.SetYmin(y_min)
		ana_unit2.SetYmax(y_max)
		ana_unit2.SetType(typ)
		ana_unit2.SetQsq(qsq)
	
		my_proc2.add_process(ana_unit2)
	
		print
		print  "Finished configuring ana_processor. Start event loop!"
		print
	
		# Let's run it.
		my_proc2.run()
	
		# done!
		print
		print "Finished running ana_processor event loop!"
		print
	
sys.exit(0)

