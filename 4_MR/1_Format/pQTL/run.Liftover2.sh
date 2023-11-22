#!/bin/bash

dir_out='aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/liftover'

arkinglab/software/src/liftOver/liftOver $dir_out/Pietzner.p05.input.round2.bed arkinglab/software/src/liftOver/hg19ToHg38.over.chain.gz $dir_out/Pietzner.p05.Liftover.Output.round2.bed $dir_out/Pietzner.p05.LiftOver.Unlifted.round2.bed
