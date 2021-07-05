#!/bin/bash

PROJ_DIR=.
DATA_DIR=${PROJ_DIR}/itc_data
OUT_DIR=${PROJ_DIR}/fits

FIT_SCRIPT=${PROJ_DIR}/../fit_itc_model.py

# Flags to run sections of the script
APPROX_FIT=1
WT_PQ1=1
C17U_PQ1=1
C31U_PQ1=1
HIN_WT_PQ1=1
NGO_WT_PQ1=1
PLOTS=0

# Number of bootstrapping iterations
N_BOOT=10000

# Plot width in inches
WIDTH=4.5

# Fit approximate model with independent sites used in PEAQ software
if [ $APPROX_FIT == 1 ]; then

DATA_1=${DATA_DIR}/can_wt_preq1_1_peaq_data
DATA_2=${DATA_DIR}/can_wt_preq1_2_peaq_data
DATA_3=${DATA_DIR}/can_wt_preq1_1_peaq_data
FIT_1=${OUT_DIR}/can_wt_preq1_1_fit_approx
FIT_2=${OUT_DIR}/can_wt_preq1_2_fit_approx
FIT_3=${OUT_DIR}/can_wt_preq1_1_fit_approx
FIT_ALL=${OUT_DIR}/can_wt_preq1_123_fit_approx
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 $DATA_1 > $FIT_1
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 $DATA_2 > $FIT_2
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 55.0 -r 3.0 -v 1420.6 $DATA_3 > $FIT_3
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 60.0,60.0,55.0 -r 3.0 -v 1420.6 $DATA_1 \
    $DATA_2 $DATA_3 > $FIT_ALL

TEMP=298.15
DATA_1=${PROJ_DIR}/can_c17u_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/can_c17u_preq1_1_peaq_data
FIT_1=${OUT_DIR}/can_c17u_preq1_1_fit_approx
FIT_2=${OUT_DIR}/can_c17u_preq1_1_fit_approx
FIT_ALL=${OUT_DIR}/can_c17u_preq1_12_fit_approx
$FIT_SCRIPT -b $N_BOOT -a -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP $DATA_1 \
    > $FIT_1
$FIT_SCRIPT -b $N_BOOT -a -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP $DATA_2 \
    > $FIT_2
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP $DATA_1 \
    $DATA_2 > $FIT_ALL

TEMP=298.15
DATA_1=${PROJ_DIR}/can_c31u_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/can_c31u_preq1_2_peaq_data
FIT_1=${OUT_DIR}/can_c31u_preq1_1_fit_approx
FIT_2=${OUT_DIR}/can_c31u_preq1_2_fit_approx
FIT_ALL=${OUT_DIR}/can_c31u_preq1_12_fit_approx
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 160.0 -r 5.0 -v 1420.6 -t $TEMP $DATA_1 \
    > $FIT_1
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 172.0 -r 4.0 -v 1420.6 -t $TEMP $DATA_2 \
    > $FIT_2
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 160.0,172.0 -r 5.0,4.0 -v 1420.6 \
    -t $TEMP $DATA_1 $DATA_2 > $FIT_ALL

DATA_1=${PROJ_DIR}/hin_wt_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/hin_wt_preq1_2_peaq_data
FIT_1=${OUT_DIR}/hin_wt_preq1_1_fit_approx
FIT_2=${OUT_DIR}/hin_wt_preq1_2_fit_approx
FIT_ALL=${OUT_DIR}/hin_wt_preq1_12_fit_approx
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 42.0 -r 3.0 -v 1420.6 $DATA_1 > $FIT_1
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 49.4 -r 3.0 -v 1420.6 $DATA_2 > $FIT_2
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 42.0,49.4 -r 3.0 -v 1420.6 $DATA_1 \
    $DATA_2 > $FIT_ALL

DATA_1=${PROJ_DIR}/ngo_wt_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/ngo_wt_preq1_2_peaq_data
FIT_1=${OUT_DIR}/ngo_wt_preq1_1_fit_approx
FIT_2=${OUT_DIR}/ngo_wt_preq1_2_fit_approx
FIT_ALL=${OUT_DIR}/ngo_wt_preq1_12_fit_approx
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 45.0 -r 3.0 -v 1420.6 $DATA_1 > $FIT_1
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 34.4 -r 2.87 -v 1420.6 $DATA_2 > $FIT_2
$FIT_SCRIPT -a -b $N_BOOT -n 2 -s 1 -l 45.0,34.4 -r 3.0,2.87 -v 1420.6 $DATA_1 \
    $DATA_2 > $FIT_ALL

fi

if [ $WT_PQ1 == 1 ]; then

E_PRIOR=16
E_STR=16
D_PRIOR=$(python3 -c 'import numpy; print(numpy.log(10))')
D_STR=log10

REG_TEST_ERROR=${OUT_DIR}/can_wt_preq1_${D_STR}_${E_STR}_reg_test_error

DATA_1=${PROJ_DIR}/can_wt_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/can_wt_preq1_2_peaq_data
DATA_3=${PROJ_DIR}/can_wt_preq1_3_peaq_data

# No regularization
PARAM_1=${OUT_DIR}/can_wt_preq1_1_no_reg_params
PARAM_2=${OUT_DIR}/can_wt_preq1_2_no_reg_params
PARAM_3=${OUT_DIR}/can_wt_preq1_3_no_reg_params

${FIT_SCRIPT} -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 $DATA_1 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_1

${FIT_SCRIPT} -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 $DATA_2 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_2

${FIT_SCRIPT} -n 2 -s 1 -l 55.0 -r 3.0 -v 1420.6 -p 0 $DATA_3 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_3

printf "no_reg %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n" \
    $($FIT_SCRIPT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_1 \
        --print_cost $DATA_2) \
    $($FIT_SCRIPT -n 2 -s 1 -l 55.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_1 \
        --print_cost $DATA_3) \
    $($FIT_SCRIPT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_2 \
        --print_cost $DATA_1) \
    $($FIT_SCRIPT -n 2 -s 1 -l 55.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_2 \
        --print_cost $DATA_3) \
    $($FIT_SCRIPT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_3 \
        --print_cost $DATA_1) \
    $($FIT_SCRIPT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_3 \
        --print_cost $DATA_2) \
    > $REG_TEST_ERROR

# Scan log_10 lambda from -6 to +6
for LOG_LAMBDA in $(seq -6 0.125 6); do

    STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
    LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

    PARAM_1=${OUT_DIR}/can_wt_preq1_1_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params
    PARAM_2=${OUT_DIR}/can_wt_preq1_2_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params
    PARAM_3=${OUT_DIR}/can_wt_preq1_3_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 60.0 -r 3.0 -v 1420.6 \
        -p $LAMBDA $DATA_1 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_1

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 60.0 -r 3.0 -v 1420.6 \
        -p $LAMBDA $DATA_2 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_2

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 55.0 -r 3.0 -v 1420.6 \
        -p $LAMBDA $DATA_3 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_3

    printf "%6.3f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n" $LOG_LAMBDA \
        $($FIT_SCRIPT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_1 \
            --print_cost $DATA_2) \
        $($FIT_SCRIPT -n 2 -s 1 -l 55.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_1 \
            --print_cost $DATA_3) \
        $($FIT_SCRIPT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_2 \
            --print_cost $DATA_1) \
        $($FIT_SCRIPT -n 2 -s 1 -l 55.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_2 \
            --print_cost $DATA_3) \
        $($FIT_SCRIPT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_3 \
            --print_cost $DATA_1) \
        $($FIT_SCRIPT -n 2 -s 1 -l 60.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_3 \
            --print_cost $DATA_2) \
        >> $REG_TEST_ERROR

done

# Find regularization penalty with lowest test error
LOG_LAMBDA=$(awk '
    NR == 1 {
        for (i = 2; i <= NF; ++i) min += $i
    } NR > 1 {
        sum = 0
        for (i = 2; i <= NF; ++i) sum += $i
        if (sum < min) {
            min = sum
            lambda = $1
        }
    } END {
        print lambda
    }' $REG_TEST_ERROR)

STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

FIT_1=${OUT_DIR}/can_wt_preq1_1_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_2=${OUT_DIR}/can_wt_preq1_2_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_3=${OUT_DIR}/can_wt_preq1_3_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_ALL=${OUT_DIR}/can_wt_preq1_123_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}

# Fit with lowest test error regularization penalty
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 60.0 -r 3.0 \
    -v 1420.6 -p $LAMBDA $DATA_1 > $FIT_1
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 60.0 -r 3.0 \
    -v 1420.6 -p $LAMBDA $DATA_2 > $FIT_2
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 55.0 -r 3.0 \
    -v 1420.6 -p $LAMBDA $DATA_3 > $FIT_3
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 60.0,60.0,55.0 \
    -r 3.0 -v 1420.6 -p $LAMBDA $DATA_1 $DATA_2 $DATA_3 > $FIT_ALL

fi

if [ $C17U_PQ1 == 1 ]; then

E_PRIOR=1
E_STR=1
D_PRIOR=1
D_STR=1
TEMP=298.15

REG_TEST_ERROR=${OUT_DIR}/can_c17u_preq1_${D_STR}_${E_STR}_reg_test_error

DATA_1=${PROJ_DIR}/can_c17u_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/can_c17u_preq1_2_peaq_data

# No regularization
PARAM_1=${OUT_DIR}/can_c17u_preq1_1_no_reg_params
PARAM_2=${OUT_DIR}/can_c17u_preq1_2_no_reg_params

${FIT_SCRIPT} -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP -p 0 $DATA_1 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_1

${FIT_SCRIPT} -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP -p 0 $DATA_2 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_2

printf "no_reg %11.6f %11.6f\n" \
    $($FIT_SCRIPT -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP -p 0 -g $PARAM_1 \
        --print_cost $DATA_2) \
    $($FIT_SCRIPT -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP -p 0 -g $PARAM_2 \
        --print_cost $DATA_1) \
    > $REG_TEST_ERROR

# Scan log_10 lambda from -6 to +6
for LOG_LAMBDA in $(seq -6 0.125 6); do

    STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
    LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

    PARAM_1=${OUT_DIR}/can_c17u_preq1_1_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params
    PARAM_2=${OUT_DIR}/can_c17u_preq1_2_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 81.0 -r 3.0 -v 1420.6 \
        -t $TEMP -p $LAMBDA $DATA_1 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_1

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 81.0 -r 3.0 -v 1420.6 \
        -t $TEMP -p $LAMBDA $DATA_2 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_2

    printf "%6.3f %11.6f %11.6f\n" $LOG_LAMBDA \
        $($FIT_SCRIPT -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP -p 0 \
            -g $PARAM_1 --print_cost $DATA_2) \
        $($FIT_SCRIPT -n 2 -s 1 -l 81.0 -r 3.0 -v 1420.6 -t $TEMP -p 0 \
            -g $PARAM_1 --print_cost $DATA_2) \
        >> $REG_TEST_ERROR

done

# Find regularization penalty with lowest test error
LOG_LAMBDA=$(awk '
    NR == 1 {
        for (i = 2; i <= NF; ++i) min += $i
    } NR > 1 {
        sum = 0
        for (i = 2; i <= NF; ++i) sum += $i
        if (sum < min) {
            min = sum
            lambda = $1
        }
    } END {
        print lambda
    }' $REG_TEST_ERROR)

STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

FIT_1=${OUT_DIR}/can_c17u_preq1_1_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_2=${OUT_DIR}/can_c17u_preq1_2_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_ALL=${OUT_DIR}/can_c17u_preq1_12_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}

# Fit with lowest test error regularization penalty
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 81.0 -r 3.0 \
    -v 1420.6 -t $TEMP -p $LAMBDA $DATA_1 > $FIT_1
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 81.0 -r 3.0 \
    -v 1420.6 -t $TEMP -p $LAMBDA $DATA_2 > $FIT_2
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 81.0 -r 3.0 \
    -v 1420.6 -t $TEMP -p $LAMBDA $DATA_1 $DATA_2 > $FIT_ALL

fi

if [ $C31U_PQ1 == 1 ]; then

E_PRIOR=16
E_STR=16
D_PRIOR=$(python3 -c 'import numpy; print(numpy.log(10))')
D_STR=log10
TEMP=298.15

REG_TEST_ERROR=${OUT_DIR}/can_c31u_preq1_${D_STR}_${E_STR}_reg_test_error

DATA_1=${PROJ_DIR}/can_c31u_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/can_c31u_preq1_2_peaq_data

# No regularization
PARAM_1=${OUT_DIR}/can_c31u_preq1_1_no_reg_params
PARAM_2=${OUT_DIR}/can_c31u_preq1_2_no_reg_params

${FIT_SCRIPT} -n 2 -s 1 -l 160.0 -r 5.0 -v 1420.6 -t $TEMP -p 0 $DATA_1 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_1

${FIT_SCRIPT} -n 2 -s 1 -l 172.0 -r 4.0 -v 1420.6 -t $TEMP -p 0 $DATA_2 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_2

printf "no_reg %11.6f %11.6f\n" \
    $($FIT_SCRIPT -n 2 -s 1 -l 172.0 -r 4.0 -v 1420.6 -t $TEMP -p 0 -g $PARAM_1 \
        --print_cost $DATA_2) \
    $($FIT_SCRIPT -n 2 -s 1 -l 160.0 -r 5.0 -v 1420.6 -t $TEMP -p 0 -g $PARAM_2 \
        --print_cost $DATA_1) \
    > $REG_TEST_ERROR

# Scan log_10 lambda from -6 to +6
for LOG_LAMBDA in $(seq -6 0.125 6); do

    STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
    LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

    PARAM_1=${OUT_DIR}/can_c31u_preq1_1_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params
    PARAM_2=${OUT_DIR}/can_c31u_preq1_2_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 160.0 -r 5.0 -v 1420.6 \
        -t $TEMP -p $LAMBDA $DATA_1 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_1

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 172.0 -r 4.0 -v 1420.6 \
        -t $TEMP -p $LAMBDA $DATA_2 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_2

    printf "%6.3f %11.6f %11.6f\n" $LOG_LAMBDA \
        $($FIT_SCRIPT -n 2 -s 1 -l 172.0 -r 4.0 -v 1420.6 -t $TEMP -p 0 \
            -g $PARAM_1 --print_cost $DATA_2) \
        $($FIT_SCRIPT -n 2 -s 1 -l 160.0 -r 5.0 -v 1420.6 -t $TEMP -p 0 \
            -g $PARAM_2 --print_cost $DATA_1) \
        >> $REG_TEST_ERROR

done

# Find regularization penalty with lowest test error
LOG_LAMBDA=$(awk '
    NR == 1 {
        for (i = 2; i <= NF; ++i) min += $i
    } NR > 1 {
        sum = 0
        for (i = 2; i <= NF; ++i) sum += $i
        if (sum < min) {
            min = sum
            lambda = $1
        }
    } END {
        print lambda
    }' $REG_TEST_ERROR)

STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

FIT_1=${OUT_DIR}/can_c31u_preq1_1_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_2=${OUT_DIR}/can_c31u_preq1_2_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_ALL=${OUT_DIR}/can_c31u_preq1_12_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}

# Fit with lowest test error regularization penalty
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 160.0 -r 5.0 \
    -v 1420.6 -t $TEMP -p $LAMBDA $DATA_1 > $FIT_1
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 172.0 -r 4.0 \
    -v 1420.6 -t $TEMP -p $LAMBDA $DATA_2 > $FIT_2
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 160.0,172.0 \
    -r 5.0,4.0 -v 1420.6 -t $TEMP -p $LAMBDA $DATA_1 $DATA_2 > $FIT_ALL

fi

if [ $HIN_WT_PQ1 == 1 ]; then

E_PRIOR=16
E_STR=16
D_PRIOR=$(python3 -c 'import numpy; print(numpy.log(10))')
D_STR=log10

REG_TEST_ERROR=${OUT_DIR}/hin_wt_preq1_${D_STR}_${E_STR}_reg_test_error

DATA_1=${PROJ_DIR}/hin_wt_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/hin_wt_preq1_2_peaq_data

# No regularization
PARAM_1=${OUT_DIR}/hin_wt_preq1_1_no_reg_params
PARAM_2=${OUT_DIR}/hin_wt_preq1_2_no_reg_params

${FIT_SCRIPT} -n 2 -s 1 -l 42.0 -r 3.0 -v 1420.6 -p 0 $DATA_1 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_1

${FIT_SCRIPT} -n 2 -s 1 -l 49.4 -r 3.0 -v 1420.6 -p 0 $DATA_2 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_2

printf "no_reg %11.6f %11.6f\n" \
    $($FIT_SCRIPT -n 2 -s 1 -l 49.4 -r 3.0 -v 1420.6 -p 0 -g $PARAM_1 \
        --print_cost $DATA_2) \
    $($FIT_SCRIPT -n 2 -s 1 -l 42.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_2 \
        --print_cost $DATA_1) \
    > $REG_TEST_ERROR

# Scan log_10 lambda from -6 to +6
for LOG_LAMBDA in $(seq -6 0.125 6); do

    STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
    LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

    PARAM_1=${OUT_DIR}/hin_wt_preq1_1_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params
    PARAM_2=${OUT_DIR}/hin_wt_preq1_2_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 42.0 -r 3.0 -v 1420.6 \
        -p $LAMBDA $DATA_1 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_1

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 49.4 -r 3.0 -v 1420.6 \
        -p $LAMBDA $DATA_2 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_2

    printf "%6.3f %11.6f %11.6f\n" $LOG_LAMBDA \
        $($FIT_SCRIPT -n 2 -s 1 -l 49.4 -r 3.0 -v 1420.6 -p 0 -g $PARAM_1 \
            --print_cost $DATA_2) \
        $($FIT_SCRIPT -n 2 -s 1 -l 42.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_2 \
            --print_cost $DATA_1) \
        >> $REG_TEST_ERROR

done

# Find regularization penalty with lowest test error
LOG_LAMBDA=$(awk '
    NR == 1 {
        for (i = 2; i <= NF; ++i) min += $i
    } NR > 1 {
        sum = 0
        for (i = 2; i <= NF; ++i) sum += $i
        if (sum < min) {
            min = sum
            lambda = $1
        }
    } END {
        print lambda
    }' $REG_TEST_ERROR)

STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

FIT_1=${OUT_DIR}/hin_wt_preq1_1_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_2=${OUT_DIR}/hin_wt_preq1_2_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_ALL=${OUT_DIR}/hin_wt_preq1_12_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}

# Fit with lowest test error regularization penalty
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 42.0 -r 3.0 \
    -v 1420.6 -p $LAMBDA $DATA_1 > $FIT_1
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 49.4 -r 3.0 \
    -v 1420.6 -p $LAMBDA $DATA_2 > $FIT_2
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 42.0,49.4 -r 3.0 \
    -v 1420.6 -p $LAMBDA $DATA_1 $DATA_2 > $FIT_ALL

fi

if [ $NGO_WT_PQ1 == 1 ]; then

E_PRIOR=1
E_STR=1
D_PRIOR=1
D_STR=1

REG_TEST_ERROR=${OUT_DIR}/ngo_wt_preq1_${D_STR}_${E_STR}_reg_test_error

DATA_1=${PROJ_DIR}/ngo_wt_preq1_1_peaq_data
DATA_2=${PROJ_DIR}/ngo_wt_preq1_2_peaq_data

# No regularization
PARAM_1=${OUT_DIR}/ngo_wt_preq1_1_no_reg_params
PARAM_2=${OUT_DIR}/ngo_wt_preq1_2_no_reg_params

${FIT_SCRIPT} -n 2 -s 1 -l 45.0 -r 3.0 -v 1420.6 -p 0 $DATA_1 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_1

${FIT_SCRIPT} -n 2 -s 1 -l 34.4 -r 2.87 -v 1420.6 -p 0 $DATA_2 | awk '
    (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
        print $3}' > $PARAM_2

printf "no_reg %11.6f %11.6f\n" \
    $($FIT_SCRIPT -n 2 -s 1 -l 34.4 -r 2.87 -v 1420.6 -p 0 -g $PARAM_1 \
        --print_cost $DATA_2) \
    $($FIT_SCRIPT -n 2 -s 1 -l 45.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_2 \
        --print_cost $DATA_1) \
    > $REG_TEST_ERROR

# Scan log_10 lambda from -6 to +6
for LOG_LAMBDA in $(seq -6 0.125 6); do

    STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
    LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

    PARAM_1=${OUT_DIR}/ngo_wt_preq1_1_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params
    PARAM_2=${OUT_DIR}/ngo_wt_preq1_2_${D_STR}_${E_STR}_reg_${STR_LOG_LAMBDA}_params

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 45.0 -r 3.0 -v 1420.6 \
        -p $LAMBDA $DATA_1 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_1

    ${FIT_SCRIPT} -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 34.4 -r 2.87 -v 1420.6 \
        -p $LAMBDA $DATA_2 | awk '
        (NR >= 5 && NR <= 8) || NR == 10 || NR == 15 || NR == 17 || NR == 18 {
            print $3}' > $PARAM_2

    printf "%6.3f %11.6f %11.6f\n" $LOG_LAMBDA \
        $($FIT_SCRIPT -n 2 -s 1 -l 34.4 -r 2.87 -v 1420.6 -p 0 -g $PARAM_1 \
            --print_cost $DATA_2) \
        $($FIT_SCRIPT -n 2 -s 1 -l 45.0 -r 3.0 -v 1420.6 -p 0 -g $PARAM_2 \
            --print_cost $DATA_1) \
        >> $REG_TEST_ERROR

done

# Find regularization penalty with lowest test error
LOG_LAMBDA=$(awk '
    NR == 1 {
        for (i = 2; i <= NF; ++i) min += $i
    } NR > 1 {
        sum = 0
        for (i = 2; i <= NF; ++i) sum += $i
        if (sum < min) {
            min = sum
            lambda = $1
        }
    } END {
        print lambda
    }' $REG_TEST_ERROR)

STR_LOG_LAMBDA=$(echo $LOG_LAMBDA | sed 's/-/n/')
LAMBDA=$(python3 -c "import numpy; print(numpy.power(10, $LOG_LAMBDA))")

FIT_1=${OUT_DIR}/ngo_wt_preq1_1_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_2=${OUT_DIR}/ngo_wt_preq1_2_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}
FIT_ALL=${OUT_DIR}/ngo_wt_preq1_12_${D_STR}_${E_STR}_fit_${STR_LOG_LAMBDA}

# Fit with lowest test error regularization penalty
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 45.0 -r 3.0 \
    -v 1420.6 -p $LAMBDA $DATA_1 > $FIT_1
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 34.4 -r 2.87 \
    -v 1420.6 -p $LAMBDA $DATA_2 > $FIT_2
$FIT_SCRIPT -b $N_BOOT -n 2 -s 1 -d $D_PRIOR -e $E_PRIOR -l 45.0,34.4 \
    -r 3.0,2.87 -v 1420.6 -p $LAMBDA $DATA_1 $DATA_2 > $FIT_ALL

fi

if [ $PLOTS == 1 ]; then

for SYSTEM in can_wt_preq1_1 can_wt_preq1_2 can_wt_preq1_3 \
    can_c17u_preq1_1 can_c17u_preq1_2 can_c31u_preq1_1 can_c31u_preq1_2 \ hin_wt_preq1_1 hin_wt_preq1_2 ngo_wt_preq1_1 ngo_wt_preq1_2; do

    FIT=${OUT_DIR}/${SYSTEM}_fit_approx

    python3 << EOF
import matplotlib.pyplot as pyplot
import numpy

pyplot.style.use('dark_background')

(molar_ratio, data, fit, confidence_lower, confidence_upper, bootstrap,
    bootstrap_lower, bootstrap_upper) = numpy.loadtxt("$FIT", unpack = True)

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio, confidence_lower, confidence_upper,
    alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio, fit, color = 'tab:orange', label = 'Fit ind sites')
pyplot.plot(molar_ratio, data, color = 'tab:blue', label = 'Data',
    linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}.pdf')

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio, bootstrap_lower, bootstrap_upper, alpha = 0.5,
    color = '0.8')
pyplot.plot(molar_ratio, bootstrap, color = 'tab:orange',
    label = 'Fit ind sites')
pyplot.plot(molar_ratio, data, color = 'tab:blue', label = 'Data',
    linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}_bootstrap.pdf')
EOF

done

for SYSTEM in can_wt_preq1_1_log10_16_fit_0.875 \
    can_wt_preq1_2_log10_16_fit_0.875 can_wt_preq1_3_log10_16_fit_0.875 \
    can_c17u_preq1_1_1_1_fit_n3.375 can_c17u_preq1_2_1_1_fit_n3.375 \
    can_c31u_preq1_1_log10_16_fit_0.375 can_c31u_preq1_2_log10_16_fit_0.375 \
    hin_wt_preq1_1_log10_16_fit_1.500 hin_wt_preq1_2_log10_16_fit_1.500 \
    ngo_wt_preq1_1_1_1_fit_n2.625 ngo_wt_preq1_2_1_1_fit_n2.625; do

    FIT=${OUT_DIR}/${SYSTEM}

    python3 << EOF
import matplotlib.pyplot as pyplot
import numpy

pyplot.style.use('dark_background')

(molar_ratio, data, fit, confidence_lower, confidence_upper, bootstrap,
    bootstrap_lower, bootstrap_upper) = numpy.loadtxt("$FIT", unpack = True)

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio, confidence_lower, confidence_upper,
    alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio, fit, color = 'tab:orange', label = 'Fit two sites reg')
pyplot.plot(molar_ratio, data, color = 'tab:blue', label = 'Data',
    linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}.pdf')

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio, bootstrap_lower, bootstrap_upper, alpha = 0.5,
    color = '0.8')
pyplot.plot(molar_ratio, bootstrap, color = 'tab:orange', \
    label = 'Fit two sites reg')
pyplot.plot(molar_ratio, data, color = 'tab:blue', label = 'Data',
    linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}_bootstrap.pdf')
EOF

done

for SYSTEM in can_c17u_preq1_12 can_c31u_preq1_12 hin_wt_preq1_12 \
    ngo_wt_preq1_12; do

    N1=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{n = split($1, s, "_")
        print (n > 3 ? s[1] "_" s[2] : s[1]) "_" s[n - 1] "_" substr(s[n], 1, 1)
        }')_fit_approx)
    N2=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{n = split($1, s, "_")
        print (n > 3 ? s[1] "_" s[2] : s[1]) "_" s[n - 1] "_" substr(s[n], 2, 1)
        }')_fit_approx)

    FIT=${OUT_DIR}/${SYSTEM}_fit_approx

    python3 << EOF
import matplotlib.pyplot as pyplot
import numpy

pyplot.style.use('dark_background')

(molar_ratio, data, fit, confidence_lower, confidence_upper, bootstrap,
    bootstrap_lower, bootstrap_upper) = numpy.loadtxt("$FIT", unpack = True)

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio[:$N1], confidence_lower[:$N1],
    confidence_upper[:$N1], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[$N1:], confidence_lower[$N1:],
    confidence_upper[$N1:], alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio[:$N1], fit[:$N1], color = 'tab:orange',
    label = 'Fit ind sites')
pyplot.plot(molar_ratio[$N1:], fit[$N1:], color = 'tab:orange', label = '')
pyplot.plot(molar_ratio[:$N1], data[:$N1], color = 'tab:blue',
    label = 'Data exp 1', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[$N1:], data[$N1:], color = 'tab:green',
    label = 'Data exp 2', linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}.pdf')

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio[:$N1], bootstrap_lower[:$N1],
    bootstrap_upper[:$N1], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[$N1:], bootstrap_lower[$N1:],
    bootstrap_upper[$N1:], alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio[:$N1], bootstrap[:$N1], color = 'tab:orange',
    label = 'Fit ind sites')
pyplot.plot(molar_ratio[$N1:], bootstrap[$N1:], color = 'tab:orange',
    label = '')
pyplot.plot(molar_ratio[:$N1], data[:$N1], color = 'tab:blue',
    label = 'Data exp 1', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[$N1:], data[$N1:], color = 'tab:green',
    label = 'Data exp 2', linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}_bootstrap.pdf')
EOF

done

for SYSTEM in can_c17u_preq1_12_1_1_fit_n3.375 \
    can_c31u_preq1_12_log10_16_fit_0.375 hin_wt_preq1_12_log10_16_fit_1.500 \
    ngo_wt_preq1_12_1_1_fit_n2.625; do

    N1=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{n = split($1, s, "_")
        print (n > 7 ? s[1] "_" s[2] : s[1]) "_" s[n - 5] "_" substr(s[n - 4], 1, 1)
        }')_fit_approx)
    N2=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{n = split($1, s, "_")
        print (n > 7 ? s[1] "_" s[2] : s[1]) "_" s[n - 5] "_" substr(s[n - 4], 2, 1)
        }')_fit_approx)

    FIT=${OUT_DIR}/${SYSTEM}

    python3 << EOF
import matplotlib.pyplot as pyplot
import numpy

pyplot.style.use('dark_background')

(molar_ratio, data, fit, confidence_lower, confidence_upper, bootstrap,
    bootstrap_lower, bootstrap_upper) = numpy.loadtxt("$FIT", unpack = True)

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio[:$N1], confidence_lower[:$N1],
    confidence_upper[:$N1], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[$N1:], confidence_lower[$N1:],
    confidence_upper[$N1:], alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio[:$N1], fit[:$N1], color = 'tab:orange',
    label = 'Fit two sites reg')
pyplot.plot(molar_ratio[$N1:], fit[$N1:], color = 'tab:orange', label = '')
pyplot.plot(molar_ratio[:$N1], data[:$N1], color = 'tab:blue',
    label = 'Data exp 1', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[$N1:], data[$N1:], color = 'tab:green',
    label = 'Data exp 2', linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}.pdf')

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio[:$N1], bootstrap_lower[:$N1],
    bootstrap_upper[:$N1], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[$N1:], bootstrap_lower[$N1:],
    bootstrap_upper[$N1:], alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio[:$N1], bootstrap[:$N1], color = 'tab:orange',
    label = 'Fit two sites reg')
pyplot.plot(molar_ratio[$N1:], bootstrap[$N1:], color = 'tab:orange',
    label = '')
pyplot.plot(molar_ratio[:$N1], data[:$N1], color = 'tab:blue',
    label = 'Data exp 1', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[$N1:], data[$N1:], color = 'tab:green',
    label = 'Data exp 2', linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}_bootstrap.pdf')
EOF

done

for SYSTEM in can_wt_preq1_123; do

    N1=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{split($1, s, "_"); print s[1] "_" s[2] "_" substr(s[3], 1, 1)
        }')_fit_approx)
    N2=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{split($1, s, "_"); print s[1] "_" s[2] "_" substr(s[3], 2, 1)
        }')_fit_approx)
    N3=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{split($1, s, "_"); print s[1] "_" s[2] "_" substr(s[3], 3, 1)
        }')_fit_approx)

    FIT=${OUT_DIR}/${SYSTEM}_fit_approx

    python3 << EOF
import matplotlib.pyplot as pyplot
import numpy

pyplot.style.use('dark_background')

i1 = $N1
i2 = $N1 + $N2

(molar_ratio, data, fit, confidence_lower, confidence_upper, bootstrap,
    bootstrap_lower, bootstrap_upper) = numpy.loadtxt("$FIT", unpack = True)

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio[:i1], confidence_lower[:i1],
    confidence_upper[:i1], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[i1:i2], confidence_lower[i1:i2],
    confidence_upper[i1:i2], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[i2:], confidence_lower[i2:],
    confidence_upper[i2:], alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio[:i1], fit[:i1], color = 'tab:orange',
    label = 'Fit ind sites')
pyplot.plot(molar_ratio[i1:i2], fit[i1:i2], color = 'tab:orange',
    label = '')
pyplot.plot(molar_ratio[i2:], fit[i2:], color = 'tab:orange', label = '')
pyplot.plot(molar_ratio[:i1], data[:i1], color = 'tab:blue',
    label = 'Data exp 1', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[i1:i2], data[i1:i2], color = 'tab:green',
    label = 'Data exp 2', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[i2:], data[i2:], color = 'tab:pink',
    label = 'Data exp 3', linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}.pdf')

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio[:i1], bootstrap_lower[:i1],
    bootstrap_upper[:i1], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[i1:i2], bootstrap_lower[i1:i2],
    bootstrap_upper[i1:i2], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[i2:], bootstrap_lower[i2:],
    bootstrap_upper[i2:], alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio[:i1], bootstrap[:i1], color = 'tab:orange',
    label = 'Fit ind sites')
pyplot.plot(molar_ratio[i1:i2], bootstrap[i1:i2], color = 'tab:orange',
    label = '')
pyplot.plot(molar_ratio[i2:], bootstrap[i2:], color = 'tab:orange',
    label = '')
pyplot.plot(molar_ratio[:i1], data[:i1], color = 'tab:blue',
    label = 'Data exp 1', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[i1:i2], data[i1:i2], color = 'tab:green',
    label = 'Data exp 2', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[i2:], data[i2:], color = 'tab:pink',
    label = 'Data exp 3', linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}_bootstrap.pdf')
EOF

done

for SYSTEM in can_wt_preq1_123_log10_16_fit_0.875; do

    N1=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{split($1, s, "_"); print s[1] "_" s[2] "_" substr(s[3], 1, 1)
        }')_fit_approx)
    N2=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{split($1, s, "_"); print s[1] "_" s[2] "_" substr(s[3], 2, 1)
        }')_fit_approx)
    N3=$(awk '$1 != "#" {++n} END {print n}' ${OUT_DIR}/$(echo $SYSTEM \
        | awk '{split($1, s, "_"); print s[1] "_" s[2] "_" substr(s[3], 3, 1)
        }')_fit_approx)

    FIT=${OUT_DIR}/${SYSTEM}

    python3 << EOF
import matplotlib.pyplot as pyplot
import numpy

pyplot.style.use('dark_background')

i1 = $N1
i2 = $N1 + $N2

(molar_ratio, data, fit, confidence_lower, confidence_upper, bootstrap,
    bootstrap_lower, bootstrap_upper) = numpy.loadtxt("$FIT", unpack = True)

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio[:i1], confidence_lower[:i1],
    confidence_upper[:i1], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[i1:i2], confidence_lower[i1:i2],
    confidence_upper[i1:i2], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[i2:], confidence_lower[i2:],
    confidence_upper[i2:], alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio[:i1], fit[:i1], color = 'tab:orange',
    label = 'Fit two sites reg')
pyplot.plot(molar_ratio[i1:i2], fit[i1:i2], color = 'tab:orange',
    label = '')
pyplot.plot(molar_ratio[i2:], fit[i2:], color = 'tab:orange', label = '')
pyplot.plot(molar_ratio[:i1], data[:i1], color = 'tab:blue',
    label = 'Data exp 1', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[i1:i2], data[i1:i2], color = 'tab:green',
    label = 'Data exp 2', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[i2:], data[i2:], color = 'tab:pink',
    label = 'Data exp 3', linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}.pdf')

figure, axes = pyplot.subplots(figsize = tuple($WIDTH * x for x in (1, 0.75)))
pyplot.fill_between(molar_ratio[:i1], bootstrap_lower[:i1],
    bootstrap_upper[:i1], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[i1:i2], bootstrap_lower[i1:i2],
    bootstrap_upper[i1:i2], alpha = 0.5, color = '0.8')
pyplot.fill_between(molar_ratio[i2:], bootstrap_lower[i2:],
    bootstrap_upper[i2:], alpha = 0.5, color = '0.8')
pyplot.plot(molar_ratio[:i1], bootstrap[:i1], color = 'tab:orange',
    label = 'Fit two sites reg')
pyplot.plot(molar_ratio[i1:i2], bootstrap[i1:i2], color = 'tab:orange',
    label = '')
pyplot.plot(molar_ratio[i2:], bootstrap[i2:], color = 'tab:orange',
    label = '')
pyplot.plot(molar_ratio[:i1], data[:i1], color = 'tab:blue',
    label = 'Data exp 1', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[i1:i2], data[i1:i2], color = 'tab:green',
    label = 'Data exp 2', linestyle = '', marker = 'D')
pyplot.plot(molar_ratio[i2:], data[i2:], color = 'tab:pink',
    label = 'Data exp 3', linestyle = '', marker = 'D')
pyplot.xlabel('Molar ratio (Ligand / Receptor)')
pyplot.ylabel('Heat per mole injected (kcal mol\$^{-1}\$)')
axes.legend(bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2)
pyplot.savefig('${FIT}_bootstrap.pdf')
EOF

done

fi

