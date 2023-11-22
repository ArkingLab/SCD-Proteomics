#!/bin/bash

dir_out='aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/liftover'

arkinglab/software/src/liftOver/liftOver $dir_out/Pietzner.p05.input.bed arkinglab/software/src/liftOver/hg19ToHg38.over.chain.gz $dir_out/Pietzner.p05.Liftover.Output.bed $dir_out/Pietzner.p05.LiftOver.Unlifted.bed
