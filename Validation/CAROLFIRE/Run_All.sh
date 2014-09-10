#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_01.fds
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_02.fds
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_03.fds
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_04.fds
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_05.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_06.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_07.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_08.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_09.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_10.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_11.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_12.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_13.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_14.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_15.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_16.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_17.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_18.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_19.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_20.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_21.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_22.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_23.fds
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_24.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_25.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_26.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_27.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_28.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_29.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_30.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_31.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_62.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_63.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_64.fds 
$QFDS $DEBUG $QUEUE -d $INDIR CAROLFIRE_PT_65.fds 

echo FDS cases submitted
