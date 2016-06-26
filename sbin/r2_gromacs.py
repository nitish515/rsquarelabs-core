__author__ = 'rrmerugu'
__VERSION__ = "0.1dev"


"""
This module provides the implementation of 'gromacs' tool step by step with the help of followed by available
commands and saves the successful project details accordingly to the project directory.

Usage:
init - initiates the project.
hello - greetings.
help - provides the available commands to be executed.
importfiles - gathers the simulation files to the project directory
createtopology -
createwaterbox -
neutralisecomplex -
minimize -

Note:
Command 'init' should not be executed in the project directory
Command 'importfiles' should be executed in the project directory
"""

# from optparse import OptionParser
from termcolor import colored, cprint
import sys, os, json, requests


from datetime import datetime
from termcolor import cprint

"""
adds the rsquarelabs-core module to this script path to access the modules inside rsquarelabs-core
"""

BIN_DIR = os.path.dirname(os.path.abspath(__file__))
CORE_DIR = os.path.join(BIN_DIR, '../')
"""
Path appended of rsquarelabs_core to sys for accessing modules inside rsquarelabs_core
"""
sys.path.append(CORE_DIR)




"""
rsquarelabs_core should be imported after the CORE_DIR is added to sys.path
"""
from rsquarelabs_core.config import RSQ_PROJECTS_HOME, RSQ_DB_PATH
from rsquarelabs_core.engines.db_engine import DBEngine
from rsquarelabs_core.engines.projects import Project
from rsquarelabs_core.engines.gromacs.gromacs import ProteinLigMin, ProteinMin
from rsquarelabs_core.engines.gromacs.core.messages import  welcome_message




"""
Used to check if the command is executed in side a project or not.
If the command is executed inside a project, 'init' will be disabled and the rest will be active.
"""
CURRENT_PATH = os.getcwd()

TOOL_NAME = "r2_gromacs"
db_object = DBEngine(RSQ_DB_PATH)


def current_date():
    """
    This method returns date and time of initiated project.

    :return: returns (str) date and time.
    """
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def show_commands():
    """

    This method provides commands for processing the project using gromacs tool and prints the available commands


    """
    available_commands = ['init', 'help', 'importfiles', 'createtopology', 'createwaterbox', 'neutralisecomplex', 'minimize']
    print "Available commands : \n"
    for command in available_commands:
        print command

def main():
    # Get the arguments list
    cmdargs = sys.argv

    # Check if config file exist in the working dir.

    files_list = os.listdir(CURRENT_PATH)
    is_config_file_avaliable = False


    for file in files_list:
        if file == "r2_gromacs.config":
            is_config_file_avaliable = True
            project_key = CURRENT_PATH.split('/')[-1]
            project_id = db_object.do_select("select id from projects where slug= ?", (project_key, )).fetchone()[0]


    if not 'init' in cmdargs:


        # Creating a object to the ProteinLigMin class


        if len(cmdargs) == 1:
            print "ERROR: Follow the allowed commands"
            show_commands()
            exit()




        # TODO - read the .yaml or project data and decide which protocol type is this
        project_data = {} # some data from db

        # lets say the protMin protocol is 1
        project_data['protocol_type'] = 1

        if project_data['protocol_type'] == 1:
            obj = ProteinMin(
                working_dir="%s/" % CURRENT_PATH,
                project_id=project_id
            )

        elif project_data['protocol_type'] == 2: # protein-lig min
            obj = ProteinLigMin(
                ligand_file='ligand.gro',
                ligand_topology_file='ligand.itp',
                protein_file='protein.pdb',
                working_dir="%s/"%CURRENT_PATH,
                project_id=project_id
            )

    #
    if 'init' in cmdargs:
        if is_config_file_avaliable:
            print "ERROR! You can't start project in this directory"
            exit()

        project_data = {}
        project_data["title"] = ""
        project_data["tags"] = ""
        project_data["user_email"] = ""
        project_data["short_note"] = ""
        project_data["slug"] = ""
        project_data["path"] = ""
        project_data["protocol_type"] = ""
        project_data["type"] = TOOL_NAME
        is_delete = 0

#         while (project_data["protocol_type"].lstrip() == ""):
#             project_data["protocol_type"] = raw_input(  """ %s
# Please select the protocol you want to work on
#
# 1. Protein Minimisation
# 2. Protein-Ligand Minimisation
# """ %(welcome_message))
#


        print "Lets start the project with Protein Minimisation Protocol"
        print "(PS: Protein-Ligand Min Protocol will be released soon)"









        while( project_data["title"].lstrip() == ""):
            project_data["title"] = raw_input("What would be your project Name: (TAK1 Modelling): ")

        while(project_data["tags"].lstrip() == ""):
            project_data["tags"] = raw_input("Please tag your project (eg: Molecular Dynamics, Minimisation, TAK1, GPCR5): ")

        while(project_data["user_email"].lstrip() == ""):
            project_data["user_email"] = raw_input("Your email for notification (eg: me@university.edu ): ")

        while (project_data["short_note"].lstrip() == ""):
            project_data["short_note"] = raw_input("Write a short note : ")

        generated_project_id = project_data["title"].replace(" ","-").replace("_","-")\
            .replace("/","-").replace("\\","-").replace(".","-").replace(",","-").replace(";",'-').replace(":","-").replace("--","-")


        customize_name = raw_input("Creating this project id as [%s], Do you wish to change ? (y/n , default=n): "%generated_project_id)


        if customize_name.lower() == 'n' or customize_name == '':
            project_data["slug"] = generated_project_id

        else:
            while(project_data["slug"].lstrip() == ""):
                project_data["slug"]  = raw_input("Enter the project_key : (tak1-modelling-trail1)")
        project_data["date"] = datetime.now().strftime("%Y-%m-%d %H:%M")

        # join rsq proj home + slug
        PROJECT_PATH = os.path.join(RSQ_PROJECTS_HOME, project_data["slug"])


        if os.path.exists(PROJECT_PATH):
            while(os.path.exists(PROJECT_PATH)):
                project_data["slug"] = raw_input("Project with project key exists, Enter new key for the project : ")
                project_data["slug"] = project_data["slug"].replace(" ","-").replace("_","-")\
                .replace("/","-").replace("\\","-").replace(".","-").replace(",","-").replace(";",'-').replace(":","-").replace("--","-")
                PROJECT_PATH = os.path.join(RSQ_PROJECTS_HOME, project_data["slug"])



        # preprocessing data


        #
        #
        # cur = db_object.do_insert("INSERT INTO projects (title, tags, user_email, slug, short_note, path, config, log, type, date, is_delete)\
        #                 VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        #                  (project_data["title"],
        #                    project_data["tags"],
        #                    project_data["user_email"],
        #                    project_data["slug"],
        #                    project_data["short_note"],
        #                    project_data["path"],
        #                    project_data["config"],
        #                    project_data["log"],
        #                    project_data["type"],
        #                    project_data["date"],
        #                    is_delete, ))
        create_project = Project(project_title=project_data["title"], project_tags=project_data["tags"], project_user_email=project_data["user_email"],
                                                   project_short_note=project_data["short_note"], project_slug=project_data["slug"])
        created_project_id = create_project.save()

        # project_data["log"] = os.path.join(PROJECT_PATH, 'r2_gromacs.log')
        # project_data["config"] = os.path.join(PROJECT_PATH, 'r2_gromacs.config')



#         log_config = db_object.do_select("select log, config from projects where id = ?", (created_project_id, )).fetchone()
#         project_data["log"] = log_config[0]
#         project_data["config"] = log_config[1]
#
#
#
#         fh_log = open(project_data["log"], 'w', 0755)
#         fh_config = open(project_data["config"], 'w', 0755)
#
#         if created_project_id: # if created into db
#             from random import randint
#             project_create_details = project_data # json.loads(project_data)
#             project_create_details['project_id'] = randint(1,1000)
#             fh_log.write("# RSQUARELABS-CORE v%s \n# Written by Ravi RT Merugu \n# https://github.com/rsquarelabs/rsquarelabs-core\n\n\n"%__VERSION__)
#
#             mesg = """============================================
# Project created with id '%s',
# ============================================""" % created_project_id
#             # fh_config.write(cur.lastrowid)
#             cprint(mesg, "green")
#         else:
#             os.remove(project_data["log"])
#             os.remove(project_data["config"])
#             os.rmdir(PROJECT_PATH)
#             mesg =  "ERROR \n%s " %project_data['title']
#             cprint(mesg, 'red')



    elif 'help' in cmdargs:
        show_commands()

    elif 'importfiles' in cmdargs:
        if is_config_file_avaliable:

            obj.import_files()
        else:
            print "ERROR! This directory do not have project details"

    elif 'createtopology' in cmdargs:
        obj.create_topology()

    elif 'createwaterbox' in cmdargs:
        obj.create_water_box()

    elif 'neutralisecomplex' in cmdargs:
        obj.neutralize_system()

    elif 'minimize' in cmdargs:
        obj.minimize()

    else:
        print "ERROR: "
        show_commands()



if __name__ == '__main__':
    main()