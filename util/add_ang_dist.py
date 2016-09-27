#!/usr/bin/python

inf = open('z64.a152.orig','r')
lines = inf.readlines()

outf = open('z64.a152', 'w')

for line in lines:
    str = line.rstrip()+"     1     0     0\n"
    outf.write(str)
