__author__ = 'rrmerugu'

import os, sys, subprocess
from datetime import datetime
from time import time
import bottle as bottle2
from bottle import Bottle, request, static_file, template, redirect, error



BASE_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)),'websuite')
STATIC_DIR  = os.path.join(BASE_DIR, 'static')
HTML_DIR  = os.path.join(STATIC_DIR, 'html')
DOCS_DIR    = os.path.join(STATIC_DIR, 'docs')
CSS_DIR     = os.path.join(STATIC_DIR, 'css')
JS_DIR      = os.path.join(STATIC_DIR, 'js')

from rsquarelabs_core.utils import run_process
from rsquarelabs_core.engines.db_engine import DBEngine
from rsquarelabs_core.config import RSQ_DB_PATH, RSQ_SCRIPT_PATH, RSQ_BACKUP_PATH, RSQ_PROJECTS_HOME

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

    if "pro_id=" in qs_string:
        pro_id = qs_string.split('pro_id=')[1].split('&')[0]

    print pro_id
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

    content = open(os.path.join(HTML_DIR, 'automator.html')).read()
    return template(content, protocol_initial_data=protocol_initial_data, new_version=new_version, now=now)

@app.route('/websuite/automator.html', method='POST')
def automator_insert():
    now = datetime.now().strftime(footer_timeformat)
    name = request.forms.get("name")
    version = request.forms.get("version")
    protocol_data = request.forms.get("protocol_data")

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


    is_delete = 0

    db_object.do_insert(" INSERT INTO protocols (name, version, parent_protocol, master_protocol, protocol_data, is_delete)\
                        VALUES ('%s', '%s', '%s', '%s', '%s', '%s')" % (name, version, parent_protocol, master_protocol, protocol_data, int(is_delete)))

    protocol_id = db_object.cur.lastrowid

    file_name = os.path.join(RSQ_SCRIPT_PATH, "%s_%s_%s.py" % (protocol_id, version, time()))
    print file_name
    fh = open(file_name, 'w')
    fh.write(protocol_data)

    redirect('/websuite/automator.html')

@app.route('/websuite/protocols.html')
def protocols():
    now = datetime.now().strftime(footer_timeformat)

    qs_string = request.query_string
    protocol_id = None

    if "protocol_id=" in qs_string:
        protocol_id = qs_string.split('protocol_id=')[1].split('&')[0]


    if protocol_id != None:
        db_object.cur.execute("UPDATE protocols SET is_delete = 1 WHERE id= ?", (protocol_id,))

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



    content = open(os.path.join(HTML_DIR, 'websuite_index.html')).read()
    return template(content, all_projects_activity=all_projects_activity,
                    all_projects_count=all_projects_count,
                    recent_projects=recent_projects,
                    recent_protocols = recent_protocols,
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
        db_object.cur.execute("UPDATE projects SET is_delete = 1 WHERE id=?", (delete_id, ))
        redirect('/websuite/projects.html')

    content = open(os.path.join(HTML_DIR, 'projects.html')).read()

    return template(content, projects_list=projects_data, now=now)


@app.route('/websuite/project/:project_id')
def projects_view(project_id):
    now = datetime.now().strftime(footer_timeformat)

    project_data = db_object.do_select("SELECT  id, slug, title, short_note, tags, user_email, type, path, log, config, date from projects where id = ?", (project_id)).fetchone()
    #TODO = filter by project_id
    project_activity_data = db_object.do_select(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity where project_id = ? ORDER BY id DESC", (project_id,))


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
    return template(content, file_list=file_list_filter, project_log=project_log, project_activity_data= project_activity_data.fetchall()[:5], project_config=project_config, project_data=project_data, now=now)


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


