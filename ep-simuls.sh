#!/bin/bash

seq 1 1 114 | parallel python2 ep-main.py --w {} --t 3000 --repl 5
