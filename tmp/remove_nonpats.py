#!/usr/bin/env python
import sys
import os

# bdir='/data1/home/inukj/mysoftware/synchdp/tutorial_data'

patl=set()
for line in open('patlist.txt').readlines():
    # print(line.strip())
    patl.add(line.strip())


for x in os.walk('./'):
    pid=x[0].split('./')[1]
    # print(pid)
    if pid not in patl:
        print(pid)
        # os.remove(x[0])
        # print(pid, ' not in list')
        # os.remove(pid)