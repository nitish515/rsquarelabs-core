__author__ = 'rrmerugu'

import os, sys
from datetime import datetime
import bottle as bottle2
from bottle import Bottle, request, static_file, template, redirect



BASE_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)),'websuite')
STATIC_DIR  = os.path.join(BASE_DIR, 'static')
HTML_DIR  = os.path.join(STATIC_DIR, 'html')
DOCS_DIR    = os.path.join(STATIC_DIR, 'docs')
CSS_DIR     = os.path.join(STATIC_DIR, 'css')
JS_DIR      = os.path.join(STATIC_DIR, 'js')

from rsquarelabs_core.utils import run_process
from rsquarelabs_core.engines.db_engine import DBEngine
from rsquarelabs_core.config import RSQ_DB_PATH

db_object = DBEngine(RSQ_DB_PATH)

app = Bottle()
now = datetime.now().strftime("%Y %b, %d %H:%M:%S %p")
bottle2.TEMPLATE_PATH.insert(0,HTML_DIR)


@app.route('/websuite/rmv.html')
def rmv():
    qs_string = request.query_string
    file_name = None
    if 'file=' in qs_string:
        file_name = qs_string.split('file=')[1].split('&')[0]
        print file_name
    content = open(os.path.join(HTML_DIR, 'rmv.html')).read()
    return template(content, file_name=file_name,now=now)


@app.route('/index')
@app.route('/home')
@app.route('/')
@app.route('/websuite')
def goto_index():
    redirect('/websuite/index.html')


@app.route('/websuite/index.html')
def index():
    print HTML_DIR
    print BASE_DIR
    content = open(os.path.join(HTML_DIR, 'websuite_index.html')).read()
    return template(content, now=now)

@app.route('/docs/rmv')
def goto_rmv():
    redirect('/docs/rmv.html')

@app.route('/docs/rmv.html')
def docs():
    content =  open(os.path.join(DOCS_DIR, 'rmv.html')).read()
    return template(content,now=now)


@app.route('/websuite/projects.html')
def projects_list():
    projects_data = db_object.do_select("SELECT id, slug, title, tags, user_email, type, path, log, date from projects")
    content =  open(os.path.join(HTML_DIR, 'projects.html')).read()
    return template(content, projects_list=projects_data,now=now)


@app.route('/websuite/project/:project_id')
def projects_view(project_id):
    project_data = db_object.do_select("SELECT  id, slug, title, short_note, tags, user_email, type, path, log, config, date from projects where id = %s"%(int(project_id))).fetchone()
    #TODO = filter by project_id

    # Filter by project id.
    qs_string = request.query_string

    id_number = None

    if 'filter_id=' in qs_string:
        id_number = qs_string.split('filter_id=')[1].split('&')[0].replace("%20", ' ')

    if id_number:
        id_number += '%'
        project_activity_data = db_object.cur.execute(
            "select id, project_id, tool_name, step_no, step_name, command, pid from project_activity where project_id LIKE ? ORDER BY id DESC",
            (id_number,))
    else:
        project_activity_data = db_object.do_select(
            "select id, project_id, tool_name, step_no, step_name, command, pid from project_activity ORDER BY id DESC")


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
    project_data = db_object.do_select("SELECT  id, slug, title, short_note, tags, user_email, type, path, log, config, date from projects where id = %s" % (int(project_id))).fetchone()

    qs_string = request.query_string

    # Filter by command name and project id.
    command_name = None
    id_number = None

    if 'filter_command=' in qs_string:
        command_name = qs_string.split('filter_command=')[1].split('&')[0].replace("%20", ' ')
    elif 'filter_id=' in qs_string:
        id_number = qs_string.split('filter_id=')[1].split('&')[0].replace("%20", ' ')

    filter_commands = ['gmx pdb2gmx', 'gmx editconf', 'gmx solvate', 'gmx grompp', 'gmx genion', 'gmx mdrun', 'gmx genrestr']

    if command_name:
        command_name += '%'
        project_activity_data = db_object.cur.execute(
            "select id, project_id, tool_name, step_no, step_name, command, pid from project_activity where command LIKE ? ORDER BY id DESC",(command_name,))
    elif id_number:
        id_number += '%'
        project_activity_data = db_object.cur.execute(
            "select id, project_id, tool_name, step_no, step_name, command, pid from project_activity where project_id LIKE ? ORDER BY id DESC",
            (id_number,))
    else:
        project_activity_data = db_object.do_select(
            "select id, project_id, tool_name, step_no, step_name, command, pid from project_activity ORDER BY id DESC")

    content = open(os.path.join(HTML_DIR, 'project-status.html')).read()
    return template(content, filter_commands=filter_commands, project_activity_data=project_activity_data.fetchall(), project_data=project_data, now=now)


@app.route('/docs/filebrowser')
def goto_filebrowser():
    redirect('/docs/filebrowser.html')

@app.route('/websuite/filebrowser.html')
def filebrowser():
    content = open(os.path.join(HTML_DIR, 'file_browser.html')).read()
    folder_path = os.path.dirname(__file__)
    folder_name = os.path.split(folder_path)[1]
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    return template(content, files_data = files, folder_name = folder_name, folder_path= folder_path,now=now )


@app.error(404)
def error404(error):
    content = open(os.path.join(HTML_DIR, '404_error.html')).read()
    return template(content)

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


