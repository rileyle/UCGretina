#!/usr/bin/python

import sys

if len(sys.argv) == 2: 
    filename = sys.argv[1]
else:
    print 'Usage: add_ang_dist.py <level data filename>'

inf = open(filename,'r')
lines = inf.readlines()

outfilename = filename+'.ad.lvldata'
outf = open(outfilename, 'w')

print "Reading "+filename+", writing "+outfilename

for line in lines:
    str = line.rstrip()+"     1     0     0\n"
    outf.write(str)
