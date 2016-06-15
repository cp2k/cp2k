#!/bin/bash
wget -r -np -e robots=off --reject "index.html*" http://www.mis.mpg.de/scicomp/EXP_SUM/1_x/
