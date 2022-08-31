#!/usr/bin/env bash

srun -A s2441 --pty -t ${1:-59} bash
