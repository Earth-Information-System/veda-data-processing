#!/usr/bin/env bash

srun -A s3673 --pty --constraint="sky|cas|hasw" -t ${1:-59} bash
