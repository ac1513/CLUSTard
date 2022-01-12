#!/bin/bash

while getopts ":w:" opt; do
  case $opt in
    w)
      wf=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [[ -n $wf ]] && [[ $wf == "both" ]];
then
  echo " ";
  echo " ===========================================";
  echo "   ____ _    _   _ ____ _____             _ ";
  echo "  / ___| |  | | | / ___|_   _|_ _ _ __ __| |";
  echo " | |   | |  | | | \___ \ | |/ _\` | '__/ _\` |";
  echo " | |___| |__| |_| |___) || | (_| | | | (_| |";
  echo "  \____|_____\___/|____/ |_|\__,_|_|  \__,_|";
  echo "                                            ";
  echo " ===========================================";
  echo "  Welcome to CLUSTard, binning and analysis ";
  echo "  workflows will be run as you specified   ";
  echo "  both! 3... 2... 1...";
  echo " ";
  snakemake -n --snakefile workflow/Snakefile_bin
  snakemake -n --snakefile workflow/Snakefile_analysis
  echo " Thank you for using CLUSTard, your metagenome bins are in the directory: ";
  echo "         results/clusters/                  ";
  echo " And the analysis output is in the directory";
  echo "         analysis                 ";
  echo " Have a good day!"
elif [[ -n $wf ]] && [[ $wf == "bin" ]] || [[ $wf == "binning" ]];
then
  echo " ";
  echo " ===========================================";
  echo "   ____ _    _   _ ____ _____             _ ";
  echo "  / ___| |  | | | / ___|_   _|_ _ _ __ __| |";
  echo " | |   | |  | | | \___ \ | |/ _\` | '__/ _\` |";
  echo " | |___| |__| |_| |___) || | (_| | | | (_| |";
  echo "  \____|_____\___/|____/ |_|\__,_|_|  \__,_|";
  echo "                                            ";
  echo " ===========================================";
  echo "  Welcome to CLUSTard, the binning workflow ";
  echo "  will be run as you specified $wf          ";
  echo "  3... 2... 1...";
  echo " ";
  snakemake -n --snakefile workflow/Snakefile_bin
  echo " Thank you for using CLUSTard, your metagenome bins are in the directory: ";
  echo "         results/clusters/                  ";
  echo " Have a good day!"
elif [[ -n $wf ]] && [[ $wf == "analysis" ]] || [[ $wf == "ana" ]];
then
  echo " ";
  echo " ===========================================";
  echo "   ____ _    _   _ ____ _____             _ ";
  echo "  / ___| |  | | | / ___|_   _|_ _ _ __ __| |";
  echo " | |   | |  | | | \___ \ | |/ _\` | '__/ _\` |";
  echo " | |___| |__| |_| |___) || | (_| | | | (_| |";
  echo "  \____|_____\___/|____/ |_|\__,_|_|  \__,_|";
  echo "                                            ";
  echo " ===========================================";
  echo "  Welcome to CLUSTard, the analysis workflow ";
  echo "  will be run as you specified $wf          ";
  echo "  The location of your metagenome bins is:  ";
  echo "  TBD";
  echo "  3... 2... 1...";
  echo " ";
  snakemake -n --snakefile workflow/Snakefile_analysis
  echo " Thank you for using CLUSTard, the analysis output is in the directory:";
  echo "         analysis                  ";
  echo " Have a good day!";
else [[ -n $wf ]];
  echo " $wf is not a valid option, please specify workflow from";
  echo " 'both', 'binning' or 'analysis'  ";
fi
