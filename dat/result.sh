#!/bin/sh

rm -f res2.out
tail -1 ARDR_GA.out >>res.out
cat res.out | awk '{print $4, $5, $6}'>>res2.out
