#!/bin/bash
#gcc buckets_of_tops.cxx -lstdc++

g++ buckets_of_tops.cxx $(root-config --cflags --libs)
