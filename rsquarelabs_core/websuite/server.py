__author__ = 'rrmerugu'

import os, sys, subprocess
from datetime import datetime
from time import time
import bottle as bottle2
from bottle import Bottle, request, static_file, template, redirect, error
from yaml import dump, load


BASE_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)),'websuite')
STATIC_DIR  = os.path.join(BASE_DIR, 'static')
HTML_DIR  = os.path.join(STATIC_DIR, 'html')
DOCS_DIR    = os.path.join(STATIC_DIR, 'docs')
CSS_DIR     = os.path.join(STATIC_DIR, 'css')
JS_DIR      = os.path.join(STATIC_DIR, 'js')

from rsquarelabs_core.utils import run_process
from rsquarelabs_core.engines.db_engine import DBEngine
from rsquarelabs_core.config import RSQ_DB_PATH, RSQ_SCRIPT_PATH, RSQ_BACKUP_PATH, RSQ_PROJECTS_HOME, RSQ_EXPORT_PATH, RSQ_IMPORT_PATH




execute_template = """
import os, sys, logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# create a file handler
handler = logging.FileHandler("RSQ_DB_LOG")
handler.setLevel(logging.INFO)
# create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(handler)


#adds the rsquarelabs-core module to this script path to access the modules inside rsquarelabs-core



CORE_DIR = "/home/nitish/PycharmProjects/rsquarelabs-core"
# Path appended of rsquarelabs_core to sys for accessing modules inside rsquarelabs_core
sys.path.append(CORE_DIR)

from rsquarelabs_core.engines.gromacs.gromacs import ProteinMin

logger.info("Executing the script")

obj = CLASS_NAME(
 project_id = "PROJECT_ID",
 working_dir = "WORKING_DIR",
 receptor_file = "RECEPTOR_FILE",
 log_file = "RSQ_DB_LOG",
 protocol_id = "PROTOCOL_ID"
)


obj.import_files()

obj.write_ions_mdp()
obj.write_prod_mdp()

"""










db_object = DBEngine(RSQ_DB_PATH)

app = Bottle()

footer_timeformat = "%Y %b, %d %H:%M:%S %p"

bottle2.TEMPLATE_PATH.insert(0, HTML_DIR)


@app.error(500)
def custom500(error):
    content = open(os.path.join(HTML_DIR, '500_error.html')).read()
    now = datetime.now().strftime(footer_timeformat)
    return template(content, now=now)

@app.error(404)
def error404(error):
    now = datetime.now().strftime(footer_timeformat)
    content = open(os.path.join(HTML_DIR, '404_error.html')).read()
    return template(content)


@app.route('/websuite/automator.html')
def automator():
    now = datetime.now().strftime(footer_timeformat)
    qs_string = request.query_string
    # pro_ver = 1
    pro_id = None
    protocol_initial_data = None
    project_recreate_data = None

    if "pro_id=" in qs_string:
        pro_id = qs_string.split('pro_id=')[1].split('&')[0]


    if pro_id:
        protocol_initial_data = db_object.do_select("SELECT id, name, version, parent_protocol, master_protocol, protocol_data from protocols WHERE id =?",
                                                    (pro_id,)).fetchone()
    new_version = db_object.do_select("SELECT MAX(version) from protocols WHERE parent_protocol=?", (pro_id, )).fetchone()[0]
    print new_version
    print type(new_version)

    if new_version != None:
        new_version = int(new_version) + 1
    else:
        if pro_id != None:
            if protocol_initial_data[3] == 0:
                new_version = 2
            else:
                new_version = 1
        else:
            new_version = 1

    if "project_id=" in qs_string:
        project_recreate_data = "\n"

        project_id = qs_string.split('project_id=')[1].split('&')[0]
        projects_list = db_object.do_select("select id, title, is_delete from projects where id=?", (project_id, )).fetchall()

        max_step = \
        db_object.do_select("select MAX(step_no) from project_activity where project_id=?", (project_id,)).fetchone()[0]

        for step_no in range(1, int(max_step)+1):
            project_recreate_activity = db_object.do_select("SELECT parent_method_name, parent_method_serial from project_activity WHERE pid_status=? and step_no=? and \
                    project_id=? ORDER BY id DESC", ("done", step_no, project_id, )).fetchone()
            print project_recreate_activity
            if project_recreate_activity[1] == 1:
                project_recreate_data = project_recreate_data + str(project_recreate_activity[0]) + "\n"
            elif project_recreate_activity == None:
                project_recreate_data = None
            #print project_recreate_data

    else:
        projects_list = db_object.do_select("select id, title, is_delete from projects",()).fetchall()

    content = open(os.path.join(HTML_DIR, 'automator.html')).read()
    return template(content, protocol_initial_data=protocol_initial_data,
                    new_version=new_version,
                    projects_list=projects_list,
                    project_recreate_data=project_recreate_data,
                    now=now)

@app.route('/websuite/automator.html', method='POST')
def automator_insert():
    now = datetime.now().strftime(footer_timeformat)
    name = request.forms.get("name")
    version = request.forms.get("version")
    protocol_data = request.forms.get("protocol_data")
    project_id = request.forms.get("project_id")
    protocol_class_name = request.forms.get("protocol_class")
    receptor_file = request.forms.get("receptor_file")

    qs_string = request.query_string
    pro_id = None
    master_protocol = 0
    parent_protocol = 0

    if "pro_id=" in qs_string:
        pro_id = qs_string.split('pro_id=')[1].split('&')[0]

        master_protocol = db_object.do_select("SELECT master_protocol from protocols WHERE id=?", (pro_id,)).fetchone()[0]
        parent_protocol = pro_id

        if master_protocol == 0:
            master_protocol = pro_id

    # if "project_id=" in qs_string:
    #     project_id = qs_string.split('project_id=')[1].split('&')[0]

    is_delete = 0
    # import re
    # protocol_data = re.escape(protocol_data)


    db_object.do_insert(" INSERT INTO protocols (name, version, parent_protocol, master_protocol, protocol_data, is_delete, project_id, class_name)\
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (name, version, parent_protocol, master_protocol, protocol_data, is_delete, project_id, protocol_class_name, ))

    protocol_id = db_object.cur.lastrowid

    project_path = db_object.do_select("select path from projects where id=?", (project_id,)).fetchone()[0]
    working_dir = os.path.join(project_path, "%s_%s" % (protocol_id, version))
    os.mkdir(working_dir, 0755)

    db_object.conn.execute("UPDATE protocols SET w_dir = ? WHERE id =?",
                           (working_dir, protocol_id,))
    db_object.conn.commit()

    file_name = os.path.join(RSQ_SCRIPT_PATH, "%s_%s_%s" % (protocol_id, version, time()))
    py_file_name = "%s.py" %file_name
    log_file_name = "%s.log" % file_name

    db_object.conn.execute("UPDATE protocols SET python_file = ?, log_file = ? WHERE id =?", (py_file_name, log_file_name, protocol_id, ))
    db_object.conn.commit()

    protocol_header = execute_template.replace('RSQ_DB_LOG', log_file_name).replace('CLASS_NAME', protocol_class_name).\
        replace('WORKING_DIR', working_dir).replace('PROJECT_ID', project_id).replace('RECEPTOR_FILE', receptor_file).\
        replace('PROTOCOL_ID', str(protocol_id))


    protocol_data_new = ""

    for line in protocol_data.split("\n"):
        line = line.rstrip().lstrip()

        if line and line.endswith(")"):

            protocol_data_new = protocol_data_new + "\nobj." + line

    protocol_data = protocol_header + protocol_data_new


    fh = open(py_file_name, 'w')
    fh.write(protocol_data)

    #TODO - find better approach for %s
    redirect('/websuite/automator_run.html?protocol_id=%s&start=true' % (protocol_id))


@app.route('/websuite/automator_run.html')
def automator_run():
    now = datetime.now().strftime(footer_timeformat)

    qs_string = request.query_string
    protocol_id = None



    if "protocol_id=" in qs_string:
        protocol_id = qs_string.split('protocol_id=')[1].split('&')[0]

    if "start=" in qs_string:
        start = qs_string.split('start=')[1].split('&')[0]



    executed_protocol = db_object.do_select(
        "SELECT python_file, log_file from protocols where id=?", (protocol_id, )).fetchone()
    if start == "true":
        subprocess.Popen(['python', executed_protocol[0]], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        redirect('/websuite/automator_run.html?protocol_id=%s&start=false' % (protocol_id))
    elif start == "false":
        from time import sleep
        sleep(2)

        if os.path.isfile(executed_protocol[1]):
            fh = open(executed_protocol[1], 'r')
            executed_data = fh.read()
        else:
            executed_data = ""

    running_activity = db_object.do_select("select id from project_activity where protocol_id=? and pid_status=?", (protocol_id,"running", )).fetchall()


    if len(running_activity) > 0:
        auto_fresh = True
    else:
        auto_fresh = False

    content = open(os.path.join(HTML_DIR, 'automator_run.html')).read()
    return template(content, executed_data=executed_data, auto_refresh=auto_fresh, now=now)


@app.route('/websuite/protocols.html')
def protocols():
    now = datetime.now().strftime(footer_timeformat)

    qs_string = request.query_string
    protocol_id = None

    if "protocol_id=" in qs_string:
        protocol_id = qs_string.split('protocol_id=')[1].split('&')[0]


    if protocol_id != None:
        db_object.conn.execute("UPDATE protocols SET is_delete = 1 WHERE id= ?", (protocol_id,))
        db_object.conn.commit()

    protocol_list = db_object.do_select("SELECT id, name, version, parent_protocol , master_protocol, protocol_data, is_delete from protocols", ())

    content = open(os.path.join(HTML_DIR, 'protocols.html')).read()
    return template(content,protocol_list=protocol_list,now=now)


@app.route('/index')
@app.route('/home')
@app.route('/')
@app.route('/websuite')
def goto_index():
    now = datetime.now().strftime(footer_timeformat)
    redirect('/websuite/index.html')


@app.route('/websuite/index.html')
def index():
    now = datetime.now().strftime(footer_timeformat)

    # TODO - optimise this query - the query is to get the total number of not deleted projects
    all_projects = db_object.do_select("select id, title from projects where is_delete = 0 ",()).fetchall()
    all_projects_activity = db_object.do_select(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity  ORDER BY id DESC",()).fetchall()

    all_projects_count = len(all_projects)
    recent_projects = all_projects[:2]
    recent_protocols = db_object.do_select(
            "select id, name, version from protocols where is_delete=0 ORDER BY id DESC ",()).fetchall()[:2]

    running_activity = db_object.do_select("select  id, step_name, project_id  from project_activity where pid_status=?", ("running", )).fetchall()
    all_running_activity_count = len(running_activity)

    content = open(os.path.join(HTML_DIR, 'websuite_index.html')).read()
    return template(content, all_projects_activity=all_projects_activity,
                    all_projects_count=all_projects_count,
                    recent_projects=recent_projects,
                    recent_protocols=recent_protocols,
                    running_activity=running_activity,
                    all_running_activity_count=all_running_activity_count,
                    now=now)


@app.route('/websuite/projects.html')
def projects_list():
    now = datetime.now().strftime(footer_timeformat)
    projects_data = db_object.do_select("SELECT id, slug, title, tags, user_email, type, path, log, date, is_delete from projects where is_delete = 0", ())

    qs_string = request.query_string
    backup_id = None
    delete_id = None

    if "backup_id=" in qs_string:
        backup_id = qs_string.split('backup_id=')[1].split('&')[0]
    elif "delete_id=" in qs_string:
        delete_id = qs_string.split('delete_id=')[1].split('&')[0]

    if backup_id != None:
        project_slug = db_object.do_select("SELECT slug from projects where id=?", (backup_id, )).fetchone()[0]
        PROJECT_ORIGIN_PATH = RSQ_PROJECTS_HOME + '/' + project_slug

        subprocess.Popen(['cp', '-r', PROJECT_ORIGIN_PATH, RSQ_BACKUP_PATH], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        redirect('/websuite/projects.html')
    elif delete_id != None:
        db_object.conn.execute("UPDATE projects SET is_delete = 1 WHERE id=?", (delete_id, ))
        db_object.conn.commit()
        redirect('/websuite/projects.html')

    content = open(os.path.join(HTML_DIR, 'projects.html')).read()

    return template(content, projects_list=projects_data, now=now)

@app.route('/websuite/export.html/:project_id')
def export(project_id):
    now = datetime.now().strftime(footer_timeformat)

    content = open(os.path.join(HTML_DIR, 'export.html')).read()
    return template(content, now=now)


@app.route('/websuite/export.html/:project_id', method='POST')
def export_yaml(project_id):
    now = datetime.now().strftime(footer_timeformat)
    # qs_string = request.query_string
    #
    # if "project_id=" in qs_string:
    #     project_id = qs_string.split('project_id=')[1].split('&')[0]
    #
    # if project_id != None:
    file_mode = request.forms.get('export')
    # print file_mode
    if file_mode:
        yaml_file = os.path.join(RSQ_EXPORT_PATH, "export_%s.yaml"%(project_id))

        # export_steps_data = []
        export_project_data = db_object.do_select("select title, tags, user_email, short_note from projects where id=?", (project_id, )).fetchone()
        # export_file_data = db_object.do_select("select  file_name, file_content  from project_files where project_id=? ORDER BY id DESC", (project_id, )).fetchone()

        export_protocol_data = db_object.do_select("select id, name, class_name from protocols where project_id=?", (project_id, )).fetchall()

        protocols_data = {}
        No_protocols = len(export_protocol_data) + 1

        if file_mode == "recreate":
            No_protocols = 1

        for protocol_order in range(1, No_protocols):
            # print export_protocol_data[protocol_order][0]
            export_activity_data = db_object.do_select("select parent_method_name, parent_method_serial, command_method from project_activity where protocol_id=?", (int(export_protocol_data[protocol_order-1][0]), )).fetchall()
            export_file_data = db_object.do_select(
                "select file_name, file_content from project_files where protocol_id = ? ",
                (int(export_protocol_data[protocol_order-1][0]), )).fetchone()
            # print export_file_data

            order = 0
            activity_data = []

            for activity_order in range(1, len(export_activity_data)+1):
                if int(export_activity_data[activity_order-1][1]) == 1:
                    order += 1
                    activity_data.append({order: str(export_activity_data[activity_order-1][0])})

            protocol_data = {
                protocol_order:
                    {
                        "protocol_name": str(export_protocol_data[protocol_order-1][1]),
                        "class_name": str(export_protocol_data[protocol_order-1][2]),
                        "files": [{str(export_file_data[0]): str(export_file_data[1])}],
                        "steps": activity_data
                    }

            }

            protocols_data.update(protocol_data)



        # if file_mode == "everything":
        #     export_steps_data = db_object.do_select("select step_no, step_name, command from project_activity where project_id=?", (project_id, )).fetchall()
        # elif file_mode == "recreate":
        #     max_step = db_object.do_select("select MAX(step_no) from project_activity where project_id=?", (project_id, )).fetchone()[0]
        #     # print int(max_step)+1
        #     for step_no in range(1, int(max_step)+1):
        #         step_data = db_object.do_select("SELECT step_no, step_name, command from project_activity WHERE pid_status=? and step_no=? and \
        #                 project_id=? ORDER BY id DESC", ("done", step_no, project_id,)).fetchone()
        #         export_steps_data.append(step_data)
        #     # print export_steps_data

        data = {
            "title": str(export_project_data[0]),
            "tag": str(export_project_data[1]),
            "user_email": str(export_project_data[2]),
            "short_note": str(export_project_data[3]),
            "protocols": protocols_data
        }

        # file_data = {
        #     "files": {
        #         str(export_file_data[0]): str(export_file_data[1])
        #     }
        # }

        # data.update(file_data)

        # steps_data = {
        #         "steps":  []
        #     }
        # for step in export_steps_data:
        #     steps_data['steps'].append({int(step[0]): {str(step[1]): str(step[2])}})
        #
        # data.update(steps_data)

        with open(yaml_file, 'w') as file_handler:
            dump(data, file_handler, default_flow_style=False)

        redirect('/websuite/download/export_%s.yaml'%(project_id))


@app.route('/websuite/download/<filename:path>')
def download(filename):
    return static_file(filename, root='%s'% (RSQ_EXPORT_PATH), download=filename)

@app.route('/websuite/import.html')
def import_yaml():
    now = datetime.now().strftime(footer_timeformat)

    content = open(os.path.join(HTML_DIR, 'import.html')).read()
    return template(content, now=now)


@app.route('/websuite/import.html', method='POST')
def import_project():
    yaml_file = request.files.get("upload")

    print yaml_file.filename
    # with open(str(yaml_file.file), 'r') as file_handler:
    #     data = load(file_handler)
    yaml_file.save(RSQ_IMPORT_PATH, overwrite=True)
    import_file = os.path.join(RSQ_IMPORT_PATH, str(yaml_file.filename))

    print import_file

    with open(import_file, 'r') as file_handler:
        data = load(file_handler)
    slug = "%s_%s"%(data['title'], time())
    db_object.do_insert("INSERT INTO PROJECTS (title, tags, user_email, short_note, slug, date, is_delete)\
        VALUES (?, ?, ?, ?, ?, ?, ?) ", (data['title'], data['tags'], data['user_email'], data['short_note'], slug, datetime.now().strftime("%Y-%m-%d %H:%M"), 0))

    redirect('/websuite/import.html')


@app.route('/websuite/project/:project_id')
def projects_view(project_id):
    now = datetime.now().strftime(footer_timeformat)

    project_data = db_object.do_select("SELECT  id, slug, title, short_note, tags, user_email, type, path, log, config, date from projects where id = ?", (project_id)).fetchone()
    #TODO = filter by project_id
    project_activity_data = db_object.do_select(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity where project_id = ? ORDER BY id DESC", (project_id,))

    protocols_list = db_object.do_select("select id, name from protocols where project_id=?", (project_id, )).fetchall()

    if project_data is None:
        project_log= None
        project_config = None
        file_list_filter = None
        project_activity_data = None
    else:
        project_log = open(project_data[8], 'r').read()
        project_config = open(project_data[9], 'r').read()
        file_list = os.listdir(project_data[7])
        file_list_filter = []
        for file in file_list:
            if not file.startswith("#") and not file.endswith("#"):
                file_list_filter.append(file)
    content = open(os.path.join(HTML_DIR, 'project-view.html')).read()
    return template(content, file_list=file_list_filter, project_log=project_log, project_activity_data= project_activity_data.fetchall()[:5], project_config=project_config, project_data=project_data,protocols_list=protocols_list, now=now)


@app.route('/websuite/project/:project_id/activity')
def activity(project_id):
    now = datetime.now().strftime(footer_timeformat)

    project_data = db_object.do_select("SELECT  id, slug, title, short_note, tags, user_email, type, path, log, config, date from projects where id = ?" , (project_id)).fetchone()
    qs_string = request.query_string

    # Filter by command name and project id.
    command_name = None


    if 'filter_command=' in qs_string:
        command_name = qs_string.split('filter_command=')[1].split('&')[0].replace("%20", ' ')


    filter_commands = ['gmx pdb2gmx', 'gmx editconf', 'gmx solvate', 'gmx grompp', 'gmx genion', 'gmx mdrun', 'gmx genrestr']

    if command_name:
        command_name += '%'
        project_activity_data = db_object.cur.execute(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity where project_id= ? and command LIKE ? ORDER BY id DESC",(project_id, command_name,))

    else:
        project_activity_data = db_object.do_select(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity  where  project_id= ?  ORDER BY id DESC", (project_id))

    content = open(os.path.join(HTML_DIR, 'project-status.html')).read()
    return template(content, filter_commands=filter_commands, project_activity_data=project_activity_data.fetchall(), project_data=project_data, now=now)






filter_commands = ['gmx pdb2gmx', 'gmx editconf']


@app.route('/assets/css/<filename>')
def server_static(filename):
    return static_file(filename, root=CSS_DIR)



def server_run():
    """
    Named _run to make this component more reusable
    """
    app.run(host='localhost', port=9090, debug=False, reloader=True, liveport=9999) #, quiet=True


"""
this is depricated
"""
def server_start_cmd():
    cmd = "nohup python %s/../bin/r2_server_start.py > /dev/null & " %(os.path.dirname(os.path.dirname(__file__)) )
    run_process(cmd)


