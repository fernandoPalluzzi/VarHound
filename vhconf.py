#! /usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys

vhlaunch = os.path.join(os.getcwd(), "vhLaunch.py")
vhconf = os.path.join(os.getcwd(), "conf/paths.ini")
manifest = os.path.join(os.getcwd(), "conf/TST500C_manifest.bed")

w = open(vhlaunch + '.tmp', "w")
with open(vhlaunch, "r") as lines:
	for line in lines:
		if line.startswith("VHCONF"):
			w.write("VHCONF = " + vhconf + "\n")
		else:
			w.write(line)
w.close()
os.rename(vhlaunch + '.tmp', vhlaunch)

w = open(vhconf + '.tmp', "w")
with open(vhconf, "r") as lines:
	for line in lines:
		if line.startswith("VarHound"):
			w.write("VarHound = " + os.path.dirname(vhlaunch) + "\n")
		elif line.startswith("manifest"):
			w.write("manifest = " + manifest + "\n")
		else:
			w.write(line)
w.close()
os.rename(vhconf + '.tmp', vhconf)
