__author__ = 'rrmerugu'
import subprocess, shlex, sys, webbrowser, logging, os, shutil
from subprocess import  Popen, PIPE
from rsquarelabs_core.engines.db_engine import  DBEngine
from rsquarelabs_core.config import RSQ_DB_PATH, RSQ_LOG_PATH
from datetime import datetime
from time import sleep

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


def set_file_premissions(file_path):

    sleep(1)
    if os.path.exists(file_path):
        os.chmod(file_path , 0777)
    else:
        print "File '%s' does not exist" %file_path
        exit()


def import_files(file_path, project_path, project_id, run_id):
    """
    This will
    1. import the files into working dir
    2. enter the file info into the database


    :param file_path: path of the importing file (from path)
    :param project_path:  path of the working dir (to path)
    :param project_id: project path
    :return:
    """
    if not os.path.exists(file_path):
        logging.error("Unable to import file '%s'" % file_path)
        exit()
    file_info = get_file_info(file_path)

    db_object.do_insert("""
    INSERT INTO project_files(file_name, file_content, project_id, run_id)
    VALUES(?,?,?,?)
    """, (file_info[0], file_info[1], project_id, run_id, ))


    """ Copy the file to project path and change permissions """
    shutil.copy2(file_path, project_path)
    os.chmod(file_path, 0777)




def get_file_info(file):
    content = open(file).read()
    file_name = file.split("/")[-1]
    return [file_name, content]

def run_and_record_process(step_no, step_name, command, tool_name, log_file, project_id, run_id, parent_method_name, parent_method_serial, command_method):
    logger.info( "INFO: Attempting to execute [STEP:%s]'%s'" %(step_no, step_name))

    try:
        ## insert the command into the db with status (to_run)
        extra = None
        if ">>" in command:
            extra = " >> %s" % command.split(">>")[1].rstrip()
            command = command.split(">>")[0].rstrip()
            cmd_args = shlex.split(command)
            if "grompp" in command:
                """
                Expection for genion command of gromacs because, it needs some data from output to be used for further step.
                """
                cmd_args = shlex.split(command).append(extra)
            extra = None

        elif "<<" in command:
            """
            #ignore <<,  because this is given when user input needs to be provided
            """
            extra = " << %s"% command.split("<<")[1].rstrip()
            command = command.split("<<")[0].rstrip()
            cmd_args =  list(shlex.split(command))

        else:
            cmd_args = shlex.split(command)






        # TODO - THIS IS INSECURE VERSION , use ? way instead of %s
        cmd = 'INSERT INTO project_activity (tool_name, step_no, step_name, command, status, log_file, project_id, created_at, run_id, parent_method_name, parent_method_serial, command_method)\
         VALUES(?,?,?,?,?,?, ?, ?,?, ?, ?, ?)'



        cur = db_object.do_insert(cmd, (tool_name, step_no, step_name, command, "to_run", log_file, project_id, datetime.now(), run_id, parent_method_name, parent_method_serial, command_method))

        logger.info(log_file)
        fh_stdout = open(log_file, 'wb')
        fh_stderr = open("%s.err"%log_file, 'wb')

        # db_object.conn.execute("UPDATE project_activity SET pid_status =? where id=? ",
        #                        ("running", cur.lastrowid))
        # db_object.conn.commit()

        logger.info("------------------------------> %s" % (cur.lastrowid))

        last_run_id = cur.lastrowid

        db_object.do_update("UPDATE project_activity SET pid_status =? where id=? ",
                               ("running", last_run_id, ))

        process = subprocess.Popen(cmd_args, stdout=PIPE, stderr=PIPE, stdin=PIPE)

        if extra:
            process.stdin.write(extra)

        stdout, stderr = process.communicate()

        fh_stdout.write(stdout)
        fh_stderr.write(stderr)
        print process.pid

        logger.info("Runing the stepNo: %s, StepName: %s with process id %s"%(step_no, step_name, process.pid))
        ret = process.poll()

        print "id of data insertion ", last_run_id

        if process.pid is not None:
            # db_object.conn.execute("UPDATE project_activity SET pid =? where id=? ", (process.pid, cur.lastrowid))
            # db_object.conn.commit()
            logger.info("------------------------------> %s"%(last_run_id))
            db_object.do_update("UPDATE project_activity SET pid =? where id=? ",
                                (process.pid, last_run_id, ))

        ret_data = {}
        ret_data['pid'] = process.pid
        ret_data['stdout'] = ""
        ret_data['stderr'] = ""


        if ret == None or ret== 0:
            logger.info('Completed!')
            # db_object.conn.execute("UPDATE project_activity SET pid_status =? where pid=? ", ("done", process.pid, ))
            # db_object.conn.commit()

            db_object.do_update("UPDATE project_activity SET pid_status =? where pid=? ", ("done", process.pid, ))
            return ret_data



        else:
            print "HeadsUP: Killed by Signal"
            logger.info( "HEADS UP: Killed by signal :(" )
            # db_object.conn.execute("UPDATE project_activity SET pid_status =? where pid=? ", ("killed", process.pid,))
            # db_object.conn.commit()
            #
            db_object.do_update("UPDATE project_activity SET pid_status =? where pid=? ", ("killed", process.pid,))
            # return ret_data
            sys.exit()

    except Exception as e:
        logger.error(e)
        logger.info( "HEADS UP: Command failed")
        # db_object.conn.execute("UPDATE project_activity SET pid_status =? where pid=? ", ("failed", process.pid,))
        # db_object.conn.commit()
        #
        db_object.do_update("UPDATE project_activity SET pid_status =? where pid=? ", ("failed", process.pid,))
        sys.exit()



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
            cmd_args = list(shlex.split(command))
            cmd_args.append(extra)
        else:
            cmd_args = shlex.split(command)




        # TODO - THIS IS INSECURE VERSION , use ? way instead of %s
        cmd = 'INSERT INTO project_activity (tool_name, step_no, step_name, command, status, log_file, project_id, created_at )\
         VALUES(?,?,?,?,?,?, ?, ?)'


        cur = db_object.do_insert(cmd, (tool_name, step_no, step_name, command, "to_run", log_file, project_id, datetime.now(), ))
        logger.debug(cur)

        fh_stdout = open(log_file, 'wb')
        fh_stderr = open("%s.err"%log_file, 'wb')
        process = Popen(cmd_args, stdout=fh_stdout, stderr=fh_stderr)

        logger.info("Runing the stepNo: %s, StepName: %s with process id %s"%(step_no, step_name, process.pid))
        ret = process.poll()

        print "id of data insertion ",cur.lastrowid

        if process.pid is not None:
            db_object.conn.execute("UPDATE project_activity SET pid =?, pid_status =? where id=? ",
                                   (process.pid, "Running", cur.lastrowid))
            db_object.conn.commit()

        ret_data = {}
        ret_data['pid'] = process.pid
        ret_data['stdout'] = fh_stdout
        ret_data['stderr'] = fh_stderr

        if ret == None:
            logger.info(  'Completed!')
            db_object.conn.execute("UPDATE project_activity SET pid_status =? where pid=? ", ("Done", process.pid,))
            db_object.conn.commit()
            return ret_data


        else:
            logger.info( "HEADS UP: Killed by signal :(", -ret)
            db_object.conn.execute("UPDATE project_activity SET pid_status =? where pid=? ", ("Killed", process.pid,))
            db_object.conn.commit()
            sys.exit()

    except Exception as e:
        logger.error(e)
        logger.info( "HEADS UP: Command failed")
        db_object.conn.execute("UPDATE project_activity SET pid_status =? where pid=? ", ("Failed", process.pid,))
        db_object.conn.commit()
        sys.exit()




def browser_register():
    webbrowser.register('r2labs_browser')

def brower_open(url):
    print "Opening the link '%s' " %url
    webbrowser.open_new(url)
    return webbrowser
