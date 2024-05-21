#!/bin/bash

d=$1
cd ${d}; 
generep2 -d -y ${d}.conf; 
cd ..;
touch ${d}-done;
