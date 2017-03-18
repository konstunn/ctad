#!/bin/bash

convmv -f cp1251 -t utf8 $1 --notest

iconv -f cp1251 -t utf8 $1 > tmp.tmp
mv tmp.tmp $1

dos2unix $1
