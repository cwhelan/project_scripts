#!/bin/bash

SVDetect linking -conf svdetect.sv.conf
SVDetect filtering -conf svdetect.sv.conf
SVDetect links2SV -conf svdetect.sv.conf
SVDetect links2bed -conf svdetect.sv.conf
