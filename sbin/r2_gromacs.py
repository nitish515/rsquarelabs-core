__author__ = 'rrmerugu'


# from optparse import OptionParser

from termcolor import colored, cprint
import sys, os, json, requests
from datetime import datetime
# get argument list using sys module
sys.argv


CLIENT_KEY = "6TCYXyf4lwTh601S1NpgbhlkyYgD5OQLbUvUq9Rf"
CLIENT_SECRET = "fZZC0uZ0aaoDICMeDsA6JXcf0ztSO7HW6t3elbQ3y4MxWdM11xGEG6l2R9zRLGxjntS5NT3bG3RcHDUL0mmJoT76PLJYHFDtSrDQFw5d6zHJ5XsyaZ9kYjX84VY82CYx"

__VERSION__ = "0.1dev"



USER_HOME_FOLDER = os.getenv('HOME')
RSQ_PROJECTS_HOME = os.path.join(USER_HOME_FOLDER, 'rsquarelabsProjects')
RSQ_PROJECTS_CONFIG = os.path.join(RSQ_PROJECTS_HOME, '.config.json')
RSQ_HOME = os.path.join(USER_HOME_FOLDER, '.rsquarelabs')
RSQ_HOME_PROJECTS_LIST = os.path.join(RSQ_HOME, 'projects-list.json')

TOOL_NAME = "r2_gromacs"

def current_date():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def main():
    # Get the arguments list
    cmdargs = str(sys.argv)
    #print cmdargs

    if 'init' in cmdargs:
        print "Lets start the project"

        # check and create folder rsquarelabsProjects in $HOME
        


        current_folder =  os.getcwd()

        project_data = {}
        project_data["project_name"] = ""
        project_data["project_tags"] = ""
        project_data["project_user_email"] = ""
        project_data["project_key"] = ""
        project_data["project_path"] = os.getcwd()
        project_data["project_type"] = TOOL_NAME

        project_data['project_config'] = os.path.join(project_data["project_path"], 'r2_gromacs.json')
        project_data['project_log'] = os.path.join(project_data["project_path"], 'r2_gromacs.log')


        if os.path.exists(project_data['project_config']) :
            mesg = "ERROR: A Project already exist in this folder\n============================================="
            cprint(mesg, 'red')
            exit()
        else:
            fh = open(project_data['project_config'] , 'w', 0755)
            fh_log = open(project_data['project_log'],'w', 0755)



        while( project_data["project_name"].lstrip() == ""):
            project_data["project_name"] = raw_input("What would be your project Name: (TAK1 Modelling): ")

        while(project_data["project_tags"].lstrip() == ""):
            project_data["project_tags"] = raw_input("Please tag your project (eg: Molecular Dynamics, Minimisation, TAK1, GPCR5): ")

        while(project_data["project_user_email"].lstrip() == ""):
            project_data["project_user_email"] = raw_input("Your email for notification (eg: me@university.edu ): ")

        generated_project_id = project_data["project_name"].replace(" ","-").replace("_","-")\
            .replace("/","-").replace("\\","-").replace(".","-").replace(",","-").replace(";",'-').replace(":","-").replace("--","-")



        customize_name = raw_input("Creating this project id as [%s], Do you wish to change ? (y/n , default=n): "%generated_project_id)

        if customize_name.lower() == 'n' or customize_name == '':
            project_data["project_key"] = generated_project_id
        else:
            while(project_data["project_key"].lstrip() == ""):
                project_data["project_key"]  = raw_input("Enter the project_key : (tak1-modelling-trail1)")

        project_data["project_started"] = datetime.now().strftime("%Y-%m-%d %H:%M")




        if os.path.exists(RSQ_HOME):
            pass
        else:
            os.mkdir(RSQ_HOME, 0755)

        # now save this config info to ~/.rsquarelabs/projects.json
        if os.path.exists(RSQ_HOME_PROJECTS_LIST):
            old_data = open(RSQ_HOME_PROJECTS_LIST).read()
            projects_list_fh = open(RSQ_HOME_PROJECTS_LIST) ## dont open in write mode
        else:
            projects_list_fh = open(RSQ_HOME_PROJECTS_LIST,"w",755)
            old_data = ""





        # print old_data
        # print len(old_data)
        if len(old_data) != 0:
            projects_list_fh = open(RSQ_HOME_PROJECTS_LIST, "w", 755)
            # convert string to dict
            data_from_file = json.loads(old_data)
            old_projects_data = data_from_file['projects']
            old_projects_data.append(project_data)

            thedata = {}
            thedata['projects'] = old_projects_data
            thedata['last_update'] =  current_date()

        else:
            projects_list_fh = open(RSQ_HOME_PROJECTS_LIST,"w",755)
            thedata = {}
            thedata['projects'] = []
            thedata['last_update'] =  current_date()
            thedata['projects'].append(project_data)
            # projects_list_fh.write(json.dumps(thedata))


        # print project_data

        ## sending this info to rsquarelabs-apis
        headers = {'content-type': 'application/json'}
        req = requests.post("http://localhost:8000/restful/project/",  headers=headers, data= json.dumps(project_data))


        if req.status_code == 201:
            project_create_details = json.loads(req.text)
            projects_list_fh.write(json.dumps(thedata))
            fh_log.write("# RSQUARELABS-CORE v%s \n# Written by Ravi RT Merugu \n# https://github.com/rsquarelabs/rsquarelabs-core\n\n\n"%__VERSION__)

            mesg = """============================================
Project created with id, %s
============================================"""  %project_create_details['project_id']
            cprint(mesg, "green")
        else:
            project_create_details = json.loads(req.text)
            mesg =  "ERROR \n%s " %project_create_details
            cprint(mesg, 'red')
            os.remove(project_data['project_config'])




    else:
        print "ERROR: please try 'r2_gromacs init'"


if __name__ == '__main__':
    main()