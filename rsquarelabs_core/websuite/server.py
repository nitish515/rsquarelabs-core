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
from rsquarelabs_core.engines.projects import Project
from rsquarelabs_core.config import RSQ_DB_PATH, RSQ_SCRIPT_PATH, RSQ_BACKUP_PATH,\
    RSQ_PROJECTS_HOME, RSQ_EXPORT_PATH, RSQ_IMPORT_PATH




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
 run_id = "RUN_ID"
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
    """
    This provides to execute the runs through server.

    :return: Returns template.
    """

    now = datetime.now().strftime(footer_timeformat)

    qs_string = request.query_string
    run_id = None
    master_id = None
    protocol_initial_data = None
    project_recreate_data = None

    # Checking to "import" from the existing runs.
    if "run_id=" in qs_string:
        run_id = qs_string.split('run_id=')[1].split('&')[0]

    # Checking to "import" from the master protocols table.
    if "master_id=" in qs_string:
        master_id = qs_string.split('master_id=')[1].split('&')[0]

    # Selecting the data from the selected run to be imported.
    if run_id:
        protocol_initial_data = db_object.do_select("SELECT run_name, run_data, parent_run_id from runs WHERE run_id =?",
                                                    (run_id,)).fetchone()

    # Selecting the data from the selected protocol to be forked.
    elif master_id:
        protocol_initial_data = db_object.do_select("SELECT protocol_name, protocol_data from protocols WHERE protocol_id =?",
                                                    (master_id,)).fetchone()
    new_version = db_object.do_select("SELECT MAX(version) from runs WHERE parent_run_id=?", (run_id, )).fetchone()[0]

    # Updating the version after importing.
    if new_version != None:
        new_version = int(new_version) + 1
    else:
        if run_id != None:
            if protocol_initial_data[2] == 0:
                new_version = 2
            else:
                new_version = 1
        else:
            new_version = 1

    # Checking for the "recreate" feature to the selected project.
    if "project_id=" in qs_string:
        project_recreate_data = "\n"

        project_id = qs_string.split('project_id=')[1].split('&')[0]

        projects_list = db_object.do_select("select id, title, is_delete from projects where id=?", (project_id, )).fetchall()

        # Finding the maximum number steps executed by the selected project from project_activity table.
        max_step = \
        db_object.do_select("select MAX(step_no) from project_activity where project_id=?", (project_id,)).fetchone()[0]

        # Collecting the successful executed data of each step but in decreasing order.
        for step_no in range(1, int(max_step)+1):
            project_recreate_activity = db_object.do_select("SELECT parent_method_name, parent_method_serial from project_activity WHERE pid_status=? and step_no=? and \
                    project_id=? ORDER BY id DESC", ("done", step_no, project_id, )).fetchone()

            if project_recreate_activity[1] == 1:
                project_recreate_data = project_recreate_data + str(project_recreate_activity[0]) + "\n"
            elif project_recreate_activity == None:
                project_recreate_data = None


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
    """
    This insert run data into runs table through automator feature.
    :return: Returns template.
    """
    now = datetime.now().strftime(footer_timeformat)

    # Collecting data from the "post" method by submitting in the local server.
    run_name = request.forms.get("run_name")
    version = request.forms.get("version")
    run_data = request.forms.get("run_data")
    project_id = request.forms.get("project_id")
    protocol_class_name = request.forms.get("protocol_class")
    receptor_file = request.forms.get("receptor_file")

    qs_string = request.query_string
    run_id = None
    master_id = 0
    parent_run_id = 0

    # Noting master_id and parent_id in case of this run is imported.
    if "run_id=" in qs_string:
        run_id = qs_string.split('run_id=')[1].split('&')[0]

        master_id = db_object.do_select("SELECT master_id from runs WHERE run_id=?", (run_id,)).fetchone()[0]
        parent_run_id = run_id

        if master_id == 0:
            master_id = run_id


    is_delete = 0

    # Saving the run.
    save_run_obj = Project()
    run_id = save_run_obj.save_run(run_name=run_name, version=version, run_data=run_data,
                              parent_run_id=parent_run_id, master_id=master_id,is_delete=is_delete,
                              protocol_class_name=protocol_class_name, project_id=project_id)

    # Creating a new directory name after this run into the project directory.
    project_path = db_object.do_select("select path from projects where id=?", (project_id,)).fetchone()[0]
    working_dir = os.path.join(project_path, "%s_%s" % (run_id, version))
    os.mkdir(working_dir, 0755)

    db_object.do_update("UPDATE runs SET w_dir = ? WHERE run_id =?",
                           (working_dir, run_id,))

    # Loading ".py" contains executed data and ".log" contains logging info, into run directory.
    file_name = os.path.join(RSQ_SCRIPT_PATH, "%s_%s_%s" % (run_id, version, time()))
    py_file_name = "%s.py" %file_name
    log_file_name = "%s.log" % file_name

    db_object.do_update("UPDATE runs SET python_file = ?, log_file = ? WHERE run_id =?", (py_file_name, log_file_name, run_id, ))

    # Adding header template to the run data.
    run_header = execute_template.replace('RSQ_DB_LOG', log_file_name).replace('CLASS_NAME', protocol_class_name).\
        replace('WORKING_DIR', working_dir).replace('PROJECT_ID', project_id).replace('RECEPTOR_FILE', receptor_file).\
        replace('RUN_ID', str(run_id))

    run_data_new = ""

    for line in run_data.split("\n"):
        line = line.rstrip().lstrip()

        if line and line.endswith(")"):
            run_data_new = run_data_new + "\nobj." + line

    run_data = run_header + run_data_new


    fh = open(py_file_name, 'w')
    fh.write(run_data)

    # Redirecting to automator_run page after submitting 'execute' button.
    #TODO - find better approach for %s
    redirect('/websuite/automator_run.html?run_id=%s&start=true' % (run_id))


@app.route('/websuite/automator_run.html')
def automator_run():
    """
    This executes the automator data and shows the simulation results during runs.
    :return: Returns template
    """
    now = datetime.now().strftime(footer_timeformat)

    qs_string = request.query_string
    run_id = None

    # Checks for loading the python data (run data) and logger data into python_file and log_file respectively.
    if "run_id=" in qs_string:
        run_id = qs_string.split('run_id=')[1].split('&')[0]

    # Checks to start the simulation.
    if "start=" in qs_string:
        start = qs_string.split('start=')[1].split('&')[0]

    executed_run = db_object.do_select(
        "SELECT python_file, log_file from runs where run_id=?", (run_id, )).fetchone()

    # If start is "true", then executes the python_file and makes it "false" to free from executing again.
    if start == "true":
        subprocess.Popen(['python', executed_run[0]], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        redirect('/websuite/automator_run.html?run_id=%s&start=false' % (run_id))
    elif start == "false":
        from time import sleep
        sleep(2)

        # Reads the log_file for displaying on this page.
        if os.path.isfile(executed_run[1]):
            fh = open(executed_run[1], 'r')
            executed_data = fh.read()
        else:
            executed_data = ""

    running_activity = db_object.do_select("select id from project_activity where run_id=? and pid_status=?", (run_id,"running", )).fetchall()

    # Checks if there is running activity so that auto_fresh makes refresh the page continuously.
    # TODO - This is not working efficiently.
    if len(running_activity) > 0:
        auto_fresh = True
    else:
        auto_fresh = False

    content = open(os.path.join(HTML_DIR, 'automator_run.html')).read()
    return template(content, executed_data=executed_data, auto_refresh=auto_fresh, now=now)


@app.route('/websuite/protocols.html')
def protocols():
    """
    This feature shows the protocols detail in a table.
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)

    # Selects the whole list of protocols for listing on this page.
    protocol_list = db_object.do_select("select protocol_id, protocol_name, author, class, description, date from protocols", ()).fetchall()

    content = open(os.path.join(HTML_DIR, 'protocols.html')).read()
    return template(content,protocol_list=protocol_list,now=now)


@app.route('/websuite/protocols/new.html')
def protocols_new():
    """
    This shows the forum to insert a protocol.
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)

    # This is for just creating new route "/websuite/protocols/new.html" and reads "new.html".
    content = open(os.path.join(HTML_DIR, 'new.html')).read()
    return template(content, now=now)


@app.route('/websuite/protocols/new.html', method='POST')
def protocols_new():
    """
    This insert the new protocol into protocols table.
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)

    # Collecting data from the "post" method by submitting in the local server.
    name = request.forms.get('name')
    author = request.forms.get('author')
    description = request.forms.get('description')
    data = request.forms.get('protocol_data')

    # Insert the master protocols into "protocols" table.
    db_object.do_insert("insert into protocols (protocol_name, class, protocol_data, description, author, date)\
     values (?, ?, ?, ?, ?, ?)", (name, "ProteinMin", data, description, author, datetime.now().strftime("%Y-%m-%d %H:%M"), ))

    redirect('/websuite/protocols/new.html')

@app.route('/websuite/protocol/:protocol_id')
def protocol_view(protocol_id):
    """
    This shows the protocol data of a protocol using its protocol_id
    :param protocol_id: identification number of a protocol
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)

    # Selecting an individual protocol for displaying the data separately.
    protocol_data = db_object.do_select("select protocol_data from protocols where protocol_id=?", (protocol_id)).fetchone()[0]

    content = open(os.path.join(HTML_DIR, 'protocol-view.html')).read()
    return template(content, protocol_data=protocol_data, now=now)


@app.route('/index')
@app.route('/home')
@app.route('/')
@app.route('/websuite')
def goto_index():
    """
    This redirect to index page.
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)
    redirect('/websuite/index.html')


@app.route('/websuite/index.html')
def index():
    """
    This provides the whole status of all projects and runs of a user.
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)

    # Selecting the all projects for displaying recent projects in the index page.
    # TODO - optimise this query - the query is to get the total number of not deleted projects
    all_projects = db_object.do_select("select id, title from projects where is_delete = 0 ",()).fetchall()
    all_projects_activity = db_object.do_select(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity  ORDER BY id DESC",()).fetchall()


    all_projects_count = len(all_projects)
    recent_projects = all_projects[:2]

    # Selecting the all runs for displaying recent runs in the index page.
    recent_runs = db_object.do_select(
            "select run_id, run_name, version from runs where is_delete=0 ORDER BY run_id DESC ",()).fetchall()[:2]

    # Selecting the running activity for displaying it.
    running_activity = db_object.do_select("select  id, step_name, project_id  from project_activity where pid_status=?", ("running", )).fetchall()
    all_running_activity_count = len(running_activity)

    content = open(os.path.join(HTML_DIR, 'websuite_index.html')).read()
    return template(content, all_projects_activity=all_projects_activity,
                    all_projects_count=all_projects_count,
                    recent_projects=recent_projects,
                    recent_runs=recent_runs,
                    running_activity=running_activity,
                    all_running_activity_count=all_running_activity_count,
                    now=now)


@app.route('/websuite/projects.html')
def projects_list():
    """
    This provides the list of projects with options of delete, back up and export to others.
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)
    projects_data = db_object.do_select("SELECT id, slug, title, tags, user_email, type, path, log, date, is_delete from projects where is_delete = 0", ())

    # qs_string for activate the actions on the project with selected id.
    qs_string = request.query_string
    backup_id = None
    delete_id = None
    export_id = None

    # Checking for the action to be activated.
    if "backup_id=" in qs_string:
        backup_id = qs_string.split('backup_id=')[1].split('&')[0]
    elif "delete_id=" in qs_string:
        delete_id = qs_string.split('delete_id=')[1].split('&')[0]
    elif "export_id=" in qs_string:
        export_id = qs_string.split('export_id=')[1].split('&')[0]

    # Back-uping the project details into "RSQ_BACKUP_PATH" directory.
    if backup_id != None:
        project_slug = db_object.do_select("SELECT slug from projects where id=?", (backup_id, )).fetchone()[0]
        PROJECT_ORIGIN_PATH = RSQ_PROJECTS_HOME + '/' + project_slug

        subprocess.Popen(['cp', '-r', PROJECT_ORIGIN_PATH, RSQ_BACKUP_PATH], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        redirect('/websuite/projects.html')

    # Deleting the project, not to display it but not in the "RSQ_PROJECTS_HOME".
    elif delete_id != None:

        db_object.do_update("UPDATE projects SET is_delete = 1 WHERE id=?", (delete_id, ))
        redirect('/websuite/projects.html')

    # Exporting the project is to make download the project in yaml format.
    elif export_id !=None:
        yaml_file = os.path.join(RSQ_EXPORT_PATH, "export_%s.yaml" % (export_id))

        # Selecting the project data to export.
        export_project_data = db_object.do_select("select title, tags, user_email, short_note from projects where id=?",
                                                  (export_id,)).fetchone()

        # Selecting the run data to export.
        export_run_data = db_object.do_select("select run_id, run_name, class_name from runs where project_id=?",
                                              (export_id,)).fetchall()

        runs_data = {}
        No_runs = len(export_run_data) + 1

        # Creating yaml format for each run in the project.
        for run_order in range(1, No_runs):

            export_activity_data = db_object.do_select(
                "select parent_method_name, parent_method_serial, command_method from project_activity where run_id=?",
                (int(export_run_data[run_order - 1][0]),)).fetchall()
            export_file_data = db_object.do_select(
                "select file_name, file_content from project_files where run_id = ? ",
                (int(export_run_data[run_order - 1][0]),)).fetchone()

            order = 0
            activity_data = []

            # Ordering the activity of a run.
            for activity_order in range(1, len(export_activity_data) + 1):
                if int(export_activity_data[activity_order - 1][1]) == 1:
                    order += 1
                    activity_data.append({order: str(export_activity_data[activity_order - 1][0])})

            # Making a dictionary of run's specification.
            run_data = {
                run_order:
            {
                "run_name": str(export_run_data[run_order - 1][1]),
                "class_name": str(export_run_data[run_order - 1][2]),
                "files": [{str(export_file_data[0]): str(export_file_data[1])}],
                "steps": activity_data
            }

            }

            # Updating the run_data right after each iteration.
            runs_data.update(run_data)

        # Making a dictionary of project's specification and adding run's dictionary.
        data = {
            "title": str(export_project_data[0]),
            "tags": str(export_project_data[1]),
            "user_email": str(export_project_data[2]),
            "short_note": str(export_project_data[3]),
            "runs": runs_data
        }

        # Dumping the dictionary into yaml file.
        with open(yaml_file, 'w') as file_handler:
            dump(data, file_handler, default_flow_style=False)

        # Redirecting to make download it.
        redirect('/websuite/download/export_%s.yaml' % (export_id))

    content = open(os.path.join(HTML_DIR, 'projects.html')).read()

    return template(content, projects_list=projects_data, now=now)


@app.route('/websuite/download/<filename:path>')
def download(filename):
    """
    This permits to download the given filename in the root path.
    :param filename: filename to be download
    :return: Returns template
    """
    return static_file(filename, root='%s'% (RSQ_EXPORT_PATH), download=filename)

@app.route('/websuite/import.html')
def import_yaml():
    now = datetime.now().strftime(footer_timeformat)

    # This is to display "import.html" page.
    content = open(os.path.join(HTML_DIR, 'import.html')).read()
    return template(content, now=now)


@app.route('/websuite/import.html', method='POST')
def import_project():
    # Getting the yaml file from uploading into the server.
    yaml_file = request.files.get("upload")

    # Save the yaml file into "RSQ_IMPORT_PATH".
    yaml_file.save(RSQ_IMPORT_PATH, overwrite=True)
    import_file = os.path.join(RSQ_IMPORT_PATH, str(yaml_file.filename))

    # Collecting the data from yaml file.
    with open(import_file, 'r') as file_handler:
        data = load(file_handler)

    # Saving the imported project into "projects" table using the collected data.
    save_project = Project(project_title=data['title'], project_tags=data['tags'],
                           project_user_email=data['user_email'], project_short_note=data['short_note'])

    project_id = save_project.save()

    # Saving the runs related to imported project into "runs" and "project_files" tables from the collected data.
    for run in data['runs']:
        import_run = data['runs'][run]
        import_file = data['runs'][run]['files']
        file_name, file_content = import_file[0].items()[0]
        import_steps = data['runs'][run]['steps']
        run_data = ''
        for step in import_steps:
            step_no, step_data = step.items()[0]
            run_data = run_data + step[step_no] + '\n'
        run_id = save_project.save_run(run_name=import_run['run_name'], protocol_class_name=import_run['class_name']
                                       , run_data=run_data, project_id=project_id)
        db_object.do_insert("INSERT INTO project_files (file_name, file_content, project_id, run_id)\
                        VALUES(?, ?, ?, ?)", (file_name, file_content, project_id, run_id))

    redirect('/websuite/import.html')


@app.route('/websuite/project/:project_id')
def projects_view(project_id):
    """
    This provides a whole view of a project using its identification number along with notes. Shows the list of runs on a project
    with export and import features.
    :param project_id: identification number of a project
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)

    # Selects the project details using project id.
    project_data = db_object.do_select("SELECT  id, slug, title, short_note, tags, user_email, type, path, log, config, date from projects where id = ?", (project_id, )).fetchone()
    #TODO = filter by project_id
    project_activity_data = db_object.do_select(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity where project_id = ? ORDER BY id DESC", (project_id,)).fetchall()

    # Selects the runs list on this project.
    runs_list = db_object.do_select("select run_id, run_name, version, class_name from runs where project_id=?", (project_id, )).fetchall()

    # Selects the Notes on this project.
    run_notes_list = db_object.do_select("select run_id, note from notes", ()).fetchall()

    # Makes all the variables to None if the project_data is None or else reads the data.
    if project_data is None:
        project_log = None
        project_config = None
        project_activity_data = None
    else:
        project_log = open(project_data[8], 'r').read()
        project_config = open(project_data[9], 'r').read()

    # qs_string for exporting an individual run.
    qs_string = request.query_string
    run_id = None

    # Checks for exporting individual run
    if "run_id=" in qs_string:
        run_id = qs_string.split('run_id=')[1].split('&')[0]
        print run_id

    if run_id !=None:
        yaml_file = os.path.join(RSQ_EXPORT_PATH, "export_%s_%s.yaml" % (project_id, run_id))

        # Selecting the project data to export.
        export_project_data = db_object.do_select("select title, tags, user_email, short_note from projects where id=?", (project_id, )).fetchone()
        # Selecting the run data to export.
        export_run_data = db_object.do_select("select run_name, class_name from runs where run_id=?", (run_id, )).fetchone()


        export_activity_data = db_object.do_select(
            "select parent_method_name, parent_method_serial, command_method from project_activity where run_id=?",
            (run_id,)).fetchall()
        export_file_data = db_object.do_select(
            "select file_name, file_content from project_files where run_id = ? ",
            (run_id,)).fetchone()


        order = 0
        activity_data = []

        # Ordering the activity of a run.
        for activity_order in range(1, len(export_activity_data) + 1):
            if int(export_activity_data[activity_order - 1][1]) == 1:
                order += 1
                activity_data.append({order: str(export_activity_data[activity_order - 1][0])})

        # Making a dictionary of run's specification.
        run_data = {
            "run_name": str(export_run_data[0]),
            "run_activity": {
                "class_name": str(export_run_data[1]),
                "files": [{str(export_file_data[0]): str(export_file_data[1])}],
                "steps": activity_data
            }
        }

        # Making a dictionary of project's specification.
        data = {
            "title": str(export_project_data[0]),
            "tags": str(export_project_data[1]),
            "user_email": str(export_project_data[2]),
            "short_note": str(export_project_data[3]),
        }

        # Updating the dictionary along with run's dictionary.
        data.update(run_data)

        # Dumping the dictionary into yaml file.
        with open(yaml_file, 'w') as file_handler:
            dump(data, file_handler, default_flow_style=False)

        # Redirecting to make download it.
        redirect('/websuite/download/export_%s_%s.yaml' % (project_id, run_id))

    content = open(os.path.join(HTML_DIR, 'project-view.html')).read()
    return template(content, run_notes_list=run_notes_list, project_log=project_log, project_activity_data= project_activity_data[:5], project_config=project_config, project_data=project_data,runs_list=runs_list, now=now)


@app.route('/websuite/runs/:run_id')
def run_veiw(run_id):
    """
    This provides view of a run using its identification number with notes and latest executed activities.
    :param run_id: identification number of a run
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)
    # Selecting the desired run's directory.
    run_dir = db_object.do_select("select w_dir from runs where run_id=?", (run_id, )).fetchone()[0]

    # Making a list of files in that directory.
    file_list = os.listdir(run_dir)
    file_list_filter = []
    for file in file_list:
        if not file.startswith("#") and not file.endswith("#"):
            file_list_filter.append(file)

    # qs_string is to allow to download a particular file.
    qs_string = request.query_string
    download_file = None

    if "download_file=" in qs_string:
        download_file = qs_string.split('download_file=')[1].split('&')[0]

    # Downloads the desired file.
    if download_file != None:
        return static_file(download_file, root='%s' %(run_dir), download=download_file)
    # Select the desired run's activity for displaying
    run_activity_data = db_object.do_select(
        "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity where run_id = ? ORDER BY id DESC",
        (run_id,)).fetchall()

    content = open(os.path.join(HTML_DIR, 'run-view.html')).read()
    return template(content, file_list=file_list_filter,
                    run_id=run_id,
                    run_activity_data=run_activity_data[:5],
                    now=now)



@app.route('/websuite/runs/:run_id', method='POST')
def run_veiw_notes(run_id):
    """
    This submits the notes regarding on the run.
    :param run_id: identification number of a run
    :return:
    """
    run_notes = request.forms.get('run_notes')
    # Submitting the notes and saves into "notes" table.
    db_object.do_insert("INSERT INTO notes (note, run_id) values (?, ?)", (run_notes, run_id, ))
    redirect('/websuite/runs/%s' %(run_id))

# @app.route('/websuite/runs/:run_id')
# def download_run_file(run_id):


@app.route('/websuite/runs/:run_id/activity')
def activity(run_id):
    """
    This shows the whole activities done with details in a particular run.
    :param run_id: identification number of a run
    :return:
    """
    now = datetime.now().strftime(footer_timeformat)

    project_id = db_object.do_select(
        "select project_id from project_activity where run_id = ? ORDER BY id DESC",
        (run_id,)).fetchone()[0]

    project_data = db_object.do_select("SELECT  id, slug, title, short_note, tags, user_email, type, path, log, config, date from projects where id = ?" , (project_id, )).fetchone()
    qs_string = request.query_string

    command_name = None

    # Checks for filtering the whole run's activity.
    if 'filter_command=' in qs_string:
        command_name = qs_string.split('filter_command=')[1].split('&')[0].replace("%20", ' ')

    # Filtering the whole activity by the list of commands provided below.
    filter_commands = ['gmx pdb2gmx', 'gmx editconf', 'gmx solvate', 'gmx grompp', 'gmx genion', 'gmx mdrun', 'gmx genrestr']

    if command_name:
        command_name += '%'
        run_activity_data = db_object.cur.execute(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity where project_id= ? and command LIKE ? ORDER BY id DESC",(run_id, command_name,))

    else:
        run_activity_data = db_object.do_select(
            "select id, tool_name, step_no, step_name, command, pid, project_id from project_activity  where  project_id= ?  ORDER BY id DESC", (run_id))

    content = open(os.path.join(HTML_DIR, 'run-status.html')).read()
    return template(content, filter_commands=filter_commands, run_activity_data=run_activity_data.fetchall(), project_data=project_data, now=now)



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


