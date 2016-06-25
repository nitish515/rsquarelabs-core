from rsquarelabs_core.config import RSQ_PROJECTS_HOME, RSQ_DB_PATH
import os
from rsquarelabs_core.engines.db_engine import DBEngine
from datetime import datetime


db_object = DBEngine(RSQ_DB_PATH)

class Project(object):

    def __init__(self):
        pass

    def create(self, *args, **kwargs):
        # check if the project-key exist in db,
        # check if the project-key exist in rsquareProjects

        project_title = kwargs.get('project_title', None)
        project_tags = kwargs.get('project_tags', None)
        project_user_email = kwargs.get('project_user_email', None)
        project_slug = kwargs.get('project_slug', None)
        project_short_note = kwargs.get('project_short_note', None)
        project_path = kwargs.get('project_path', None)
        # project_config = kwargs.get('project_config', None)
        # project_log = kwargs.get('project_log', None)

        # PROJECT_PATH = os.path.join(RSQ_PROJECTS_HOME, project_slug)

        # if os.path.exists(PROJECT_PATH):
        #     #print "ERROR: Project name already exists"
        #     os.rmdir(PROJECT_PATH)
        #     os.mkdir(PROJECT_PATH, 0755)
        # else:
        #     os.mkdir(PROJECT_PATH, 0755)

        project_type ="r2_gromacs"
        project_date = datetime.now().strftime("%Y-%m-%d %H:%M")
        is_delete = 0
        # project_path = PROJECT_PATH
        project_log = os.path.join(project_path, 'r2_gromacs.log')
        project_config = os.path.join(project_path, 'r2_gromacs.config')

        db_object.do_insert("INSERT INTO projects (title, tags, user_email, slug, short_note, path, config, log, type, date, is_delete)\
            VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                   (project_title,
                                    project_tags,
                                    project_user_email,
                                    project_slug,
                                    project_short_note,
                                    project_path,
                                    project_config,
                                    project_log,
                                    project_type,
                                    project_date,
                                    is_delete,))
        # return cur
        project_id = db_object.cur.lastrowid

        return project_id



    def create_run(self, *args, **kwargs):

        run_name = kwargs.get('run_name', None)
        version = kwargs.get('version', None)
        parent_run_id = kwargs.get('parent_run_id', None)
        master_id = kwargs.get('master_id', None)
        run_data = kwargs.get('run_data', None)
        is_delete = kwargs.get('is_delete', None)
        project_id = kwargs.get('project_id', None)
        protocol_class_name = kwargs.get('protocol_class_name', None)

        db_object.do_insert(" INSERT INTO runs (run_name, version, parent_run_id, master_id, run_data, is_delete, project_id, class_name)\
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (run_name, version, parent_run_id, master_id, run_data, is_delete, project_id, protocol_class_name, ))

        run_id = db_object.cur.lastrowid

        return run_id

