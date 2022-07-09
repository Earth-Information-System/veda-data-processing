#!/usr/bin/env bash

srun -A s2441 --pty --constraint="sky|cas|hasw" -t ${1:-59} bash
