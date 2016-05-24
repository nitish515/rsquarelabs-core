__author__ = 'rrmerugu'
import subprocess, shlex, sys, webbrowser, logging
from rsquarelabs_core.engines.db_engine import  DBEngine
from rsquarelabs_core.config import RSQ_DB_PATH, RSQ_LOG_PATH
from datetime import datetime

# logging.basicConfig(filename=RSQ_LOG_PATH, level=logging.DEBUG)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# create a file handler
handler = logging.FileHandler(RSQ_LOG_PATH)
handler.setLevel(logging.DEBUG)
# create a logging format
formatter = logging.Formatter('%(asctime)s - %(lineno)d - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(handler)



db_object = DBEngine(RSQ_DB_PATH)

def check_process():
    pass


def run_process(step_no, step_name, command, tool_name, log_file, project_id):
    """
    This method will
    step1: save the incoming command request info into db
    step2: executes the command
    step3: saves the pid of running and changes the command status to executing
    (whether executed or not is checked by other method 'check_process()' )

    :param step_no: Step number in the workflow (this helps us in tracking how many times the user failed at this step
    and also we use the last of this step no from the records as the final command that worked for the user :) )
    :param step_name: Step name in the workflow - an identifier
    :param command: the command to execute
    :param project_id: id of project (this helps for filtering in project status)
    :return:
    """
    logger.info( "INFO: Attempting to execute " + step_name + " [STEP:" + step_no + "]")
    try:
        ## insert the command into the db with status (to_run)
        extra = ""
        if ">>" in command:
            extra = " >> %s" % command.split(">>")[1].rstrip()
            command = command.split(">>")[0].rstrip()
            cmd_args = shlex.split(command)
            if "grompp" in command:
                """
                Expection for genion command of gromacs because, it needs some data from output to be used for further step.
                """
                cmd_args = shlex.split(command).append(extra)


        elif "<<" in command:
            """
            #ignore <<,  because thhis is given when user input needs to be provided
            """
            extra = " << %s"% command.split("<<")[1].rstrip()
            command = command.split("<<")[0].rstrip()
            cmd_args =  list(shlex.split(command))
            cmd_args.append(extra)
        else:
            cmd_args = shlex.split(command)




        # TODO - THIS IS INSECURE VERSION , use ? way instead of %s
        cmd = 'INSERT INTO project_activity (tool_name, step_no, step_name, command, status, log_file, project_id, created_at )\
         VALUES("%s",%s,"%s","%s","%s","%s", %s, "%s")'% (tool_name, int(step_no), step_name, command, "to_run", log_file, int(project_id), str(datetime.now()) )


        cur = db_object.do_insert(cmd)
        logger.debug(cur)

        fh_stdout = open(log_file, 'wb')
        fh_stderr = open("%s.err"%log_file, 'wb')
        process = subprocess.Popen(cmd_args, stdout=fh_stdout, stderr=fh_stderr)

        logger.info("Runing the stepNo: %s, StepName: %s with process id %s"%(step_no, step_name, process.pid))
        ret = process.poll()

        print "id of data insertion ",cur.lastrowid

        if process.pid is not None:
            db_object.cur.execute( "UPDATE project_activity SET pid ='%s' where id='%s' " %(process.pid, cur.lastrowid))


        ret_data = {}
        ret_data['pid'] = process.pid
        ret_data['stdout'] = fh_stdout
        ret_data['stderr'] = fh_stderr

        if ret == None:
            logger.info(  'Completed!')
            return ret_data


        else:
            logger.info( "HEADS UP: Killed by signal :(", -ret)
            sys.exit()

    except Exception as e:
        logger.error(e)
        logger.info( "HEADS UP: Command failed")
        sys.exit()




def browser_register():
    webbrowser.register('r2labs_browser')

def brower_open(url):
    print "Opening the link '%s' " %url
    webbrowser.open_new(url)
    return webbrowser
