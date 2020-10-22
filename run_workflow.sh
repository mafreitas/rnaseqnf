#!/usr/bin/env bash
# Project: Quantiattive Proteomics Pipeline
# Developer: Michael A. Freitas
#
# Description: Nextflow, OpenMM, Python, Proteowizard, Quantiative Proteomics
# Copyright (c) 2020, Michael A. Freitas, The Ohio State University
#
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 

set -e

function usage
{
    echo "usage: run_pipeline -d -s SOME_MORE_ARGS [-y YET_MORE_ARGS || -h]"
    echo "   ";
    echo "  -h | --help : This message";
}

function parse_args
{
  # positional args
  args=()

  # named args
  while [ "$1" != "" ]; do
    case "$1" in
      -h | --help )         usage;            exit;; # quit and show usage
      * )                   args+=("$1")             # if no match, add it to the positional args
    esac
    shift 
    # move to next kv pair
  done

  # restore positional args
  set -- "${args[@]}"
}


function run
{
  parse_args "$@"
  set -e
  export TOWER_ACCESS_TOKEN='d6055031d82dc8c6d54c03f3428158ca4f96d061'
  module load nextflow
  nextflow main.nf -resume -profile slurm -with-tower
}

run "$@";
