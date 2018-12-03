#!/bin/sh
if [ $1 == 'clean' ]
then
	cd io
	pwd
	make clean
	cd ../checkpoint
	make clean
	cd ../calc
	make clean
	cd ../
else
	cd io
	make
	cd ../checkpoint
	make
	cd ../calc
	make
	cd bin
fi
