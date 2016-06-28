#!/usr/bin/env bash
#This will clean all the data files including db, projects files.

PROJ_FOLDER=$HOME"/rsquarelabsProjects/"
RSQ_HOME_FOLDER=$HOME"/.rsquarelabs/"

echo "Sanitising the RSQ_HOME_FOLDER :" $RSQ_HOME_FOLDER
rm -Rf $RSQ_HOME_FOLDER
echo "Sanitising the PROJ_FOLDER :" $PROJ_FOLDER
rm -Rf $PROJ_FOLDER

