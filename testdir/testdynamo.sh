#/bin/sh          
# W Mi, June 18 2006
# This is a script to generate banchmark data.

#clear screen 
clear

#Check OS and Hostname pick MAG Dynamo source
OS=`/bin/uname |sed 's/[0-9 . -]//g'`
HOSTNAME=`/bin/hostname`
echo "OS is ${OS}"
echo "Hostname is ${HOSTNAME}"

#if  [ $HOSTNAME = "dynamo" ]
#then
#    MAG_DIR="/home/wei/dynamo"
#else
#   if  [ $HOSTNAME = "www" ]
#   then
#     MAG_DIR="/home/wei/magdynamo"
#   fi
#fi

#change this to your source directory
MAG_DIR="/home/wei/magdynamo"
echo "MAG Daynamo source is $MAG_DIR"

echo "----------------------------------"
echo " Main Menu"
echo "----------------------------------"
echo "[1] Test Bench0"
echo "[2] Test Bench1"
echo "[3] Generate Bench0 Data Files"
echo "[4] Generate Bench1 Data Files"
echo "[5] Exit/Stop"
echo "=================================="
echo "Please select the test case for test:[1-5]"
read MAG_TEST

#Define MAG test directory
MAGTEST_DIR="$MAG_DIR/test-dynamo"
#Run_testcase()
#{
case $MAG_TEST in     
  1)         
    echo "The test is for MAG Dynamo benchmark0"
    BENCHMARK_DIR="$MAGTEST_DIR/data-bench0-olson" 
    #Define test data directory and infile
    TESTDATA_DIR="$MAGTEST_DIR/test-data-bench0" 
    INFILE="$MAG_DIR/par.bnch0"
    RUN_MAG=$MAG_DIR/magx32s4
    ;;             
  2)             
    echo "The test is for MAG Dynamo benchmark1"
    BENCHMARK_DIR="$MAGTEST_DIR/data-bench1-olson" 
    #Define test data directory
    TESTDATA_DIR="$MAGTEST_DIR/test-data-bench1"
    INFILE="$MAG_DIR/par.bnch1"    
    RUN_MAG=$MAG_DIR/magx32s4 
    ;;
  3)         
    echo "This is going to generate MAG Dynamo Benchmark0 Data"
    BENCHMARK_DIR="$MAGTEST_DIR/data-bench0"
#    echo "The Benchmark direcotory is $BENCHMARK_DIR" 
    #Define infile
    INFILE="$MAG_DIR/par.bnch0"
    RUN_MAG=$MAG_DIR/magx32s4
    $RUN_MAG <$INFILE
    mv $MAG_DIR/*.bench0 $BENCHMARK_DIR/.
    exit 0
    ;;             
  4)             
    echo "This going to generate MAG Dynamo benchmark1 Data"
    BENCHMARK_DIR="$MAGTEST_DIR/data-bench1" 
    #Define infile
    INFILE="$MAG_DIR/par.bnch1"    
    RUN_MAG=$MAG_DIR/magx32s4
    $RUN_MAG <$INFILE
    mv $MAG_DIR/*.bench1 $BENCHMARK_DIR/.
    exit 0
    ;;                     
  5)         
    exit 0
    ;; 
  *)
    exit 0
    ;;            
esac       
#}

#Clean previous testdata files 
echo "Cleaning previous testdata"
cd $TESTDATA_DIR
rm ./*.*

#Test start time
TIME_START=`date +%m.%d.%Y.%H:%M`
echo "Test started at $TIME_START"

#Run MAG Dynamo
$RUN_MAG <$INFILE

#Test finish time
TIME_FINISH=`date +%m.%d.%Y.%H:%M`
echo "Test finished at $TIME_FINISH"


#Move all output files from MAG directory to test data directory
mv $MAG_DIR/*.bench* $TESTDATA_DIR/.

#check_result()

cd $MAGTEST_DIR/test-report
echo "clean Test Report Directory"
rm ./*.diff

cd $TESTDATA_DIR

#Do diff on all output files in test data directory
#Save result in Test Report
 for i in `find -type f -print`
 do
 #
 echo "  Checking $i"
 diff $BENCHMARK_DIR/$i $i > $MAGTEST_DIR/test-report/$i.diff 
 #
 done


#Check test result
cd $MAGTEST_DIR/test-report
for i in `find *.diff -type f -print`
do
#  echo -n "."
  echo "Checking $i ..."
  s=`sum $i|cut -d" " -f1`
      if [ "$s" -eq "0" ]; then
        echo "$s:${i}"
        echo " No Difference in Data File $i"
	     echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> sum.$TIME_FINISH
        echo "No Difference found in Data File $i" >> sum.$TIME_FINISH 
        echo "                    " >> sum.$TIME_FINISH
        cat $i >> sum.$TIME_FINISH
      else
        echo "$s:${i}" 
        echo "$i: Failed"
        echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> sum.$TIME_FINISH
        echo "Test on $i Failed" >> sum.$TIME_FINISH 
        echo "                    " >> sum.$TIME_FINISH
        cat $i >> sum.$TIME_FINISH
      fi
done
echo "Test started at $TIME_START" >> sum.$TIME_FINISH
echo "Test finished at $TIME_FINISH" >> sum.$TIME_FINISH
echo "*****************************************************"
echo " Please check summary in Report directory for details" 
echo " summary file : sum.$TIME_FINISH"
echo "Test started at $TIME_START"
echo "Test finished at $TIME_FINISH"
echo "*****************************************************"


exit
